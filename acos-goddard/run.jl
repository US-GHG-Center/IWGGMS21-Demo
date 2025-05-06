# The retrieval toolkit - every function you use with a RE.xxx is from there
# Any other function is either from a third-party module or from this repo.

using RetrievalToolbox
const RE = RetrievalToolbox


using ArgParse
using Dates
using DocStringExtensions
using HDF5
using Interpolations
using LinearAlgebra
using Logging
using LoopVectorization
using Printf
using ProgressMeter
using Statistics
using StatsBase
using Unitful


# Functions to read in data
include("io.jl")
# Function to create a state vector
include("state_vector.jl")
# Function that processes a scenen through the inversion algorithm
include("process_scene.jl")
# Forward model
include("forward_model.jl")
# Non-uniform sampling
include("NUS.jl")
# Routines to process the results
include("collect_results.jl")

logger = ConsoleLogger(stderr, Logging.Info);
global_logger(logger);


#=
    Global variables
=#

# Number of pressure levels on the retrieval grid
N_RT_lev = 20
# what number type we use (keep this Float64)
my_type = Float64

# This is where the cmdline arguments are processed
include("args.jl")


function julia_main() #::Cint

    #=
        Prepare the paths that contain the input data
    =#

    # Parse the command line arguments into a dictionary
    args = parse_commandline()

    # Parse the spectral windows
    spec_array = [parse(Int, x) for x in split(args["spec"], ",")]

    # Which spectrometers are we processing?
    if 1 in spec_array
        @info "Processing O2 A-band"
    end
    if 2 in spec_array
        @info "Processing Weak CO2 band"
    end
    if 3 in spec_array
        @info "Processing Strong CO2 band"
    end

    #=
        Read in the retrieval inputs
    =#

    @info "Reading in scene inputs ..."

    scene_inputs = Dict(
        spec => OCO_read_inputs_from_l1b_and_met(
            args["L1b"],
            args["L2Met"],
            band_idx=spec,
            sounding_id=args["sounding_id"],
            l2cpr_fname=args["L2CPr"]
            ) for spec in spec_array
        )

    @info "... done!"

    # Skip bad soundings immediately
    _input = first(values(scene_inputs))
    if _input["sounding_qual_flag"][args["sounding_id"]] != 0
        @info "Bad sounding quality. Exiting."
        exit(1)
    end

    #=
        #########################
        Create the needed objects
        #########################
    =#

    # For the OCO-type solar model
    @info "Reading in solar model from $(args["solar_model"]) .. "

    if args["solar_model"][end-2:end] == ".h5"
        # Use the ACOS solar model
        solar_models = Dict(
                spec => RE.OCOHDFSolarModel(
                # Path to the solar model file
                args["solar_model"],
                # "spec" referring to OCO-HDF solar model band index
                spec,
                spectral_unit=:Wavelength
            ) for spec in spec_array
        )

    elseif args["solar_model"][end-2:end] == ".nc"
        # Use the TSIS model
        solar_models = Dict(
            spec => RE.TSISSolarModel(
                args["solar_model"],
                spectral_unit=:Wavelength
            ) for spec in spec_array
        )

        for (spec, sm) in solar_models
            # In-place conversion of solar model units
            RE.convert_solar_model_to_photons!(sm)
        end
    end
    @info "... done!"


    # Read in spectroscopy, depending on the window configuration. We
    # keep them in a Dictionary that are accessed via a simple string
    # so that e.g. `abscos["CO2"]` gets you the CO2 absco.
    abscos = Dict{String, RE.AbstractSpectroscopy}()

    # If we use the O2 A-band, we need O2 at least
    if 1 in spec_array
        @info "Reading in O2 spectroscopy ..."
        abscos["O2"] = RE.load_ABSCO_spectroscopy(
            # Pass the path to the ABSCO file
            "./data/o2_v52.hdf",
            spectral_unit=:Wavelength
        )
        @info "... done!"

        # Apply a user-defined spectroscopy scaling factor. This can be done to the
        # entire spectroscopic table since oxygen is only present in the A-band
        # anyway.
        @info "Scaling O2 cross sections by $(args["o2_scale"])"
        abscos["O2"].cross_section[:] .*= args["o2_scale"]

    end

    # If we use either of the two CO2 bands, we need H2O and CO2
    if (2 in spec_array) | (3 in spec_array)
        @info "Reading in CO2 spectroscopy ..."
        abscos["CO2"] = RE.load_ABSCO_spectroscopy(
            # Pass the path to the ABSCO file
            "./data/co2_v52.hdf",
            spectral_unit=:Wavelength
        )
        @info "... done!"

        # Apply a user-defined scale factor for CO2 spectroscopy in the weak band
        @info "Scaling CO2 cross sections for weak CO2 by " *
            "$(args["co2_scale_weak"])"
        idx_weak = findall(abscos["CO2"].ww * abscos["CO2"].ww_unit .< 2.0u"µm")
        abscos["CO2"].cross_section[idx_weak,:,:,:] .*= args["co2_scale_weak"]


        # Apply a user-defined scale factor for CO2 spectroscopy in the strong band
        @info "Scaling CO2 cross sections for strong CO2 by " *
            "$(args["co2_scale_strong"])"
        idx_strong = findall(abscos["CO2"].ww * abscos["CO2"].ww_unit .> 2.0u"µm")
        abscos["CO2"].cross_section[idx_strong,:,:,:] .*= args["co2_scale_strong"]


        @info "Reading in H2O spectroscopy ..."
        abscos["H2O"] = RE.load_ABSCO_spectroscopy(
            # Pass the path to the ABSCO file
            "./data/h2o_v52.hdf",
            spectral_unit=:Wavelength
        )
        @info "... done!"
    end


    #=
        For every ABSCO object, we create a new gas object. Same as the ABSCOs,
        we stick them into a dictionary that is accessed by the gas name, e.g.
        `gases["O2"]` delivers the oxygen gas object which itself uses the
        `abscos["O2"]` ABSCO table.
    =#


    gas_units = Dict(
        "O2" => Unitful.NoUnits,
        "H2O" => Unitful.percent,
        "CO2" => Unitful.ppm
        )

    gases = Dict{String, RE.GasAbsorber}()

    for (gas_name, gas_absco) in abscos
        @info "Creating gas $(gas_name) with ABSCO $(gas_absco.file_name) ..."

        gases[gas_name] = RE.GasAbsorber(
            gas_name, # Name of the gas
            gas_absco, # ABSCO object
            zeros(my_type, N_RT_lev), # Placeholder for VMR level profile
            gas_units[gas_name]
        )

        @info "... done!"
    end


    #=
        Create the spectral windows!
        (they are created using a helper function will will make sure
         that the spectral grid is congruent with the ABSCO sampling)
    =#

    spectral_windows = Dict{Int, RE.SpectralWindow}()

    if 1 in spec_array

        spectral_windows[1] = RE.spectralwindow_from_ABSCO(
            "o2",
            0.7592152960425568, 0.7714610963322822, 0.760,
            #0.759203, 0.771429, 0.760,
            #0.75925, 0.771, 0.760,
            5.0e-3,
            abscos["O2"],
            Unitful.μm # or u"μm"
        )

    end

    if 2 in spec_array
        spectral_windows[2] = RE.spectralwindow_from_ABSCO(
            "weak_co2",
            #1.5955, 1.6057, 1.60,
            #1.5955, 1.608, 1.60,
            #1.597, 1.62, 1.60,
            #1.597986, 1.61778, 1.60,
            1.5979699810546206, 1.617757907862207, 1.60,
            5.0e-3,
            abscos["CO2"],
            Unitful.μm # or u"μm"
        )
    end

    if 3 in spec_array
        spectral_windows[3] = RE.spectralwindow_from_ABSCO(
            "strong_co2",
            #2.048, 2.08, 2.0,
            #2.047653993018265, 2.068, 2.0,
            #2.04773, 2.079794, 2.0,
            #2.04773, 2.078, 2.0,
            2.0476103795182246, 2.0796702096182256, 2.0,
            5.0e-3,
            abscos["CO2"],
            Unitful.μm # or u"μm"
        )
    end


    h5_oco_static = h5open("./example_data/l2_oco_static_input.h5", "r")

    #=
        Non-uniform sampling (NUS) grid
        (read the NUS knots from the OCO static file)
    =#



    nus_dict = Dict{Int, Vector{Bool}}()
    for spec in spec_array
        @info "Creating non-uniform sampling data for $(spectral_windows[spec])."
        # The NUS knots are in wavenumber
        nus_ww = h5_oco_static["Spectrum_Sampling/nonuniform_grid_$(spec)"][:]
        # .. convert to wavelength
        nus_wl = 1e4 ./ nus_ww
        # find which wavelengths in this spectral window are in common with the
        # NUS knots from the file
        wl_common = intersect(nus_wl, spectral_windows[spec].ww_grid)
        # Where are the common wavelengths found?
        nus_idx = searchsortedfirst.(Ref(spectral_windows[spec].ww_grid), wl_common);

        # Now create a new boolean mask array with the same length as the
        # high-res wavelengths
        nus = zeros(Bool, spectral_windows[spec].N_hires);
        # .. and set it to `true` for every spectral point that is to be calculated!
        nus[nus_idx] .= 1;
        # Store it in the dict.
        nus_dict[spec] = nus;
    end


    # Grab the CO2 covariance matrix from the ACOS static file
    # (and divide by (ppm^2) to get to ppm units.
    CO2_covar = h5_oco_static["Gas/CO2/covariance"][:,:] .* 1e12

    close(h5_oco_static)


    #=
        Create dispersion objects!
    =#
    my_dispersions = Dict{Int, RE.SimplePolynomialDispersion}()
    for spec in spec_array
        @info "Creating dispersion object for $(spectral_windows[spec])."

        # Note that the dispersion polynomial coefficients are present
        # in every instance of `scene_inputs`.

        my_dispersions[spec] = RE.SimplePolynomialDispersion(
            scene_inputs[spec]["dispersion"][spec],
            1:1016,
            spectral_windows[spec]
        );
    end



    #=
        Create the atmospheric elements to be used in the retrieval!
    =#

    # Start with an empty vector
    atm_elements = RE.AbstractAtmosphereElement[]

    # Add all available gases
    for (gas_name, gas) in gases
        push!(atm_elements, gas)
    end

    # Add Rayleigh scattering
    push!(atm_elements, RE.RayleighScattering())

    #=
        Add the aerosols we use for every scene (ice cloud, water cloud, statospheric sulphur)
    =#

    if args["aerosols"]

        @info "Using aerosols."

        #=
            Ice cloud
        =#
        @info "Adding ice cloud ..."
        miemom_ice = read_ACOS_aerosol(
            "example_data/l2_aerosol_combined_slim.h5",
            "ice_cloud_MODIS6_deltaM_50",
            wavelength=true
        )
        aer_ice = RE.GaussAerosol(
            "Ice",
            miemom_ice,
            true, # relative pressure?
            0.0, # height p/psurf
            0.04, # width p/psurf
            Unitful.NoUnits,
            0.755, # ref wavelength
            u"µm", # ref wavelength unit
            exp(-5.382), # Total AOD at reference point
            );

        push!(atm_elements, aer_ice)

        #=
            Water cloud
        =#
        @info "Adding water cloud ..."
        miemom_water = read_ACOS_aerosol(
            "example_data/l2_aerosol_combined_slim.h5",
            "wc_008",
            wavelength=true
        )
        aer_water = RE.GaussAerosol(
            "Water",
            miemom_water,
            true, # relative pressure?
            0.75, # height p/psurf
            0.1, # width p/psurf
            Unitful.NoUnits,
            0.755, # ref wavelength
            u"µm", # ref wavelength unit
            exp(-4.382), # Total AOD at reference point
            );

        push!(atm_elements, aer_water)

        #=
            Stratospheric aerosol
        =#
        @info "Adding stratospheric aerosol ..."
        miemom_strat = read_ACOS_aerosol(
            "example_data/l2_aerosol_combined_slim.h5",
            "strat",
            wavelength=true
            );
        aer_strat = RE.GaussAerosol(
            "ST",
            miemom_strat,
            true, # relative pressure?
            0.03, # height p/psurf
            0.04, # width p/psurf
            Unitful.NoUnits,
            0.755, # ref wavelength
            u"µm", # ref wavelength unit
            exp(-5.116), # Total AOD at reference point
        );

        push!(atm_elements, aer_strat)

        #=
            Now add the two largest contributors to tropospheric aerosols from MERRA
        =#
        # These are potentially hard-coded, but really look at "/Metadata/CompositeAerosolTypes"
        trop_types = ["DU", "SS", "BC", "OC", "SO"]

        # (add 1 to account for 0-based vs 1-based indexing..)
        trop_idx = scene_inputs[first(spec_array)]["met"][args["sounding_id"]]["aerosol_sort_met"] .+ 1
        # Find where the indices are 1 or 2 (to get the top two aerosols)
        # (or adjust to your liking, ..)
        trop_choice_idx = indexin([1,2], trop_idx)
        trop_chosen = trop_types[trop_choice_idx]

        @info "Top two aerosols are: $(trop_chosen)"

        for idx in trop_choice_idx
            @info "Reading aerosol properties for $(trop_types[idx])"

            this_trop_miemom = read_ACOS_aerosol(
                "example_data/l2_aerosol_combined_slim.h5",
                trop_types[idx],
                wavelength=true
            )

            # Extract height (fraction of psurf), width (fraction of psurf), AOD
            gauss_params = scene_inputs[first(spec_array)]["met"][args["sounding_id"]]["aerosol_gauss_params_met"][:, idx]

            # The position of the Gaussian parameters was inferred from the L2FP code..
            # (might need checking)
            this_trop_height = convert(my_type, gauss_params[2])
            this_trop_width = convert(my_type, gauss_params[3])
            # Override -- ACOS says we have to do 0.05 width
            this_trop_width = 0.05
            this_trop_aod = convert(my_type, gauss_params[4])

            @info "Creating Gauss-type aerosol for $(trop_types[idx])."
            this_trop_gauss = RE.GaussAerosol(
                trop_types[idx],
                this_trop_miemom,
                true, # relative pressure?
                this_trop_height, # height p/psurf
                this_trop_width, # width p/psurf
                Unitful.NoUnits,
                0.755, # ref wavelength
                u"µm", # ref wavelength unit
                this_trop_aod, # Total AOD at reference point
            )
            # Push this aerosol into the list of atmospheric elements

            push!(atm_elements, this_trop_gauss)

        end
    end
    # Finished populating the list of atmospheric elements!
    @info "Atmospheric elements: $(atm_elements)"



    #@info "Creating state vector ..."
    # Create a state vector with a helper function.
    state_vector = create_statevector(
        spectral_windows,
        my_dispersions,
        atm_elements
    )

    N_sv = length(state_vector)

    #=

        Create buffers.
        Buffers are pre-allocated objects and arrays which the various
        functions then use to in-place modify and place intermediate
        results into. This dramatically reduces the amount of allocations
        which would negatively impact overall performance.

    =#

    @info "Create buffers ..."

    # Will contain outputs of ISRF application. These values seem to work for
    # a normal 3-band OCO-type retrieval.
    inst_buf = RE.InstrumentBuffer(
        zeros(my_type, 5000),
        zeros(my_type, 5000),
        zeros(my_type, 1016),
    )

    #=
    Will contain ouptuts of OE calculations
    (this can be optional, maybe the chosen inversion
     method does not require an OE buffer, or maybe
     it requires something else)
    =#

    oe_buf = RE.OEBuffer(
        3*1016, N_sv, my_type
        )

    # Will contain outputs of radiance, jacobians and
    # dispersion indices. We make this a ScalarRTBuffer because
    # the instrument only measures Intensity, and not the polarization components.

    # We use ScalarRTBuffer because the instrument measures scalar radiance!
    rt_buf = RE.ScalarRTBuffer(
        Dict(spectral_windows[i] => my_dispersions[i] for i in spec_array),
        RE.ScalarRadiance(my_type, 3*1016), # Hold the radiance
        Dict(sve => RE.ScalarRadiance(my_type, 3*1016) for sve in state_vector.state_vector_elements),
        Dict(spectral_windows[i] => zeros(Int, 0) for i in spec_array), # Hold the detector indices
        u"ph/s/m^2/sr/µm", # Radiance units for the forward model, this must be in the same units as the L1B data..
    )

    # Will contain everything needed to run a retrieval

    buf = RE.EarthAtmosphereBuffer(
        state_vector, # The state vector
        [spectral_windows[i] for i in spec_array], # The spectral window (or a list of them)
        #[(:RPV, 3, 1.0, -0.1, 0.05) for x in spec_array], # Surface type for each band (same order as spectral windows)
        [(:Lambertian, 3) for x in spec_array], # Surface type for each band (same order as spectral windows)
        atm_elements, # Atmospheric elements
        Dict(spectral_windows[i] => solar_models[i] for i in spec_array), # The solar model(s),
        [:XRTM for i in spec_array], # RT model
        RE.VectorRadiance, # We use VectorRadiance here because we want to do the calculations with full Stokes!
        rt_buf, # The pre-allocated RT buffer
        inst_buf, # The pre-allocated instrument buffer
        N_RT_lev, # The number of retrieval or RT pressure levels
        72, # The number of meteorological pressure levels
        my_type # The chosen Float data type (e.g. Float16, Float32, Float64)
    )

    @info "... done!"


    # Setting the RT options

    # Grab the number of high streams from ARGS
    Nhigh = 16


    # The monochromatic calculations are doing using first-order polarized RT, and then
    # a fast two-stream multiple scattering scalar calculation.
    mo1 = Dict()
    mo2 = Dict()

    mo1["solvers"] = ["single"]
    mo1["add"] = true
    mo1["streams"] = 2
    mo1["options"] = [
        "output_at_levels",
        "calc_derivs",
        "source_solar",
        "vector",
        "psa",
        "sfi",
        ]

    if args["aerosols"]
        push!(mo1["options"], "delta_m")
        push!(mo1["options"], "n_t_tms")
    end

    mo2["solvers"] = ["two_stream"]
    mo2["add"] = true
    mo2["streams"] = 2
    mo2["options"] = [
        "output_at_levels",
        "calc_derivs",
        "source_solar",
        "psa",
        "sfi"
        ]

    if args["aerosols"]
        push!(mo2["options"], "delta_m")
        push!(mo2["options"], "n_t_tms")
    end

    for (swin, rt) in buf.rt

        empty!(rt.model_options)
        push!(rt.model_options, mo1)
        push!(rt.model_options, mo2)

    end


    #=
        High-stream options for LSI
    =#

    hmo1 = Dict()
    hmo1["solvers"] = ["single"]
    hmo1["add"] = true
    hmo1["streams"] = Nhigh
    hmo1["options"] = [
        "output_at_levels",
        "calc_derivs",
        "source_solar",
        "vector",
        "psa",
        "sfi",
        ]

    if args["aerosols"]
        push!(hmo1["options"], "delta_m")
        push!(hmo1["options"], "n_t_tms")
    end

    hmo2 = Dict()
    hmo2["solvers"] = ["two_os"]
    hmo2["add"] = true
    hmo2["streams"] = Nhigh
    hmo2["options"] = [
        "output_at_levels",
        "calc_derivs",
        "source_solar",
        "vector",
        "psa",
        "sfi",
    ]

    hmo3 = Dict()
    hmo3["solvers"] = ["eig_bvp"]
    hmo3["add"] = true
    hmo3["streams"] = Nhigh
    hmo3["options"] = [
        "output_at_levels",
        "calc_derivs",
        "source_solar",
        "psa",
        "sfi",
        ]

    if args["aerosols"]
        push!(hmo3["options"], "delta_m")
        push!(hmo3["options"], "n_t_tms")
    end

    if (args["LSI"])
        high_options = [
            hmo1,
            hmo2,
            hmo3
        ]
    else
        high_options = nothing
    end

    #=
        Process the scene
    =#

    @time solver, fm_kwargs = process_snid(
        spectral_windows,
        scene_inputs,
        state_vector,
        args["sounding_id"],
        buf,
        oe_buf;
        max_iter=args["max_iterations"],
        gamma=args["gamma"],
        dsigma_scale=args["dsigma_scale"],
        co2_cov=CO2_covar,
        high_options=high_options,
        nus_dict=nus_dict,
        )

    return buf, solver, fm_kwargs
    
end

julia_main()
