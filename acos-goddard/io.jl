"""
    Walks through an HDF5 file `h5` and reads everything into Dict `d`!
"""

function walk_h5!(h5, d; root="")
    for key in keys(h5)
        if h5[key] isa HDF5.Group
            walk_h5!(h5[key], d, root="$(root)/$(key)")
        else
            d["$(root)/$(key)"] = HDF5.read_dataset(h5, key)
        end
    end
end




"""
Reads necessary retrieval inputs from NASA ACOS L1bSc, L2Met and L2CPr (optional) files.

## Details

This function reads in the whole-orbit scene and instrument data needed to run XGAS
retrievals. Many quantities read from the corresponding files are supplied with physical
units via the Unitful.jl package. The units of the quantities are not explicity stated in
the source HDF5 files, however the can (mostly) be inferred or looked up in the official
OCO documents (see below).

The results are generally stored in dictionaries, whose keys are either sounding IDs
themselves, or, when appropriate, the keys would be e.g. the OCO footprint number. For
example, the meteorological variable are part of the "met" dictionary, which itself
consists of dictionaries that contiain the needed data. To access, e.g., the surface
pressure for sounding ID 12345, one would do

oco_inputs = OCO_read_inputs_from_l1b_and_met(l1b_fname, met_fname, l2cpr_fname,
band_idx=1) this_psurf = oco_inputs["met"]["surface_pressure_met"][12345]

Note that `this_psurf` now has units of `hPa`, as per the `OCO_pull_atmosphere_from_met`
function.

OCO documents: https://disc.gsfc.nasa.gov/information/documents?title=OCO-2%20Documents

"""
function OCO_read_inputs_from_l1b_and_met(
    l1b_fname::String,
    met_fname::String;
    l2cpr_fname = nothing,
    sounding_id = nothing,
    band_idx = 1
)

    # Open L1B and L2Met files
    l1b_h5 = h5open(l1b_fname, "r", swmr=true)
    met_h5 = h5open(met_fname, "r", swmr=true)

    # Make the L2CPr optional
    if !isnothing(l2cpr_fname)
        cpr_h5 = h5open(l2cpr_fname, "r", swmr=true)
    else
        cpr_h5 = nothing
    end

    # Obtain the sounding ID list
    if isnothing(sounding_id)
        sounding_id_list = vec(l1b_h5["/SoundingGeometry/sounding_id"][:,:])
    else
        sounding_id_list = [sounding_id]
    end

    # Get the number of footprints and frames in this granule
    N_fp, N_fr = size(l1b_h5["/SoundingGeometry/sounding_id"])

    # This here is a dictionary that maps sounding id -> (footprint, frame)
    fp_frame_idx_dict = RE.OCO_get_footprint_and_frame_from_snid(l1b_h5)

    # Grab the MET data from the files
    met_dicts = RE.OCO_pull_atmosphere_from_met(met_h5)
    # Get the scene locations as location objects
    loc_dict = RE.OCO_pull_sounding_location(l1b_h5)
    # Get the sounding times
    time_dict = RE.OCO_pull_sounding_time(l1b_h5)
    # Create and get the observers (satellite angles, velocity, position)
    observer_dict = RE.OCO_pull_observer(l1b_h5)

    if !isnothing(cpr_h5)
        @info "Pulling CO2 priors from $(l2cpr_fname)"
        co2_prior_dict = RE.OCO_pull_CO2_prior(cpr_h5)
    end

    #=
        Additional scene-dependent quantities needed for the retrieval
    =#
    meas_full_dict = Dict() # sounding ID -> full-range measured radiance
    stokes_coefs_dict = Dict() # sounding ID -> Stokes coefficients per scene
    sza_dict = Dict() # sounding ID -> SZA per scene
    saa_dict = Dict() # sounding ID -> SAA per scene
    sounding_qual_dict = Dict() # sounding ID -> sounding quality flag
    doppler_factor_dict = Dict() # sounding ID -> instrument Doppler factor (or better: term)
    solar_doppler_factor_dict = Dict() # sounding ID -> solar Doppler factor (or better: term)
    solar_distance_dict = Dict() # sounding ID -> normalized solar distance

    #=
        Instrument parameters, these are the same for scene in a given footprint
    =#
    isrf_table_dict = Dict() # footprint -> ISRF table
    disp_coef_dict = Dict() # footprint -> dispersion coeffs
    snr_coef_dict = Dict() # footprint -> per-wl SNR coefficient arrays

    @info "Reading ISRF and dispersion tables"
    for fp in 1:N_fp

        isrf_table_dict[fp] = RE.TableISRF(
            l1b_h5["InstrumentHeader/ils_delta_lambda"][:,:,fp,band_idx],
            u"µm",
            l1b_h5["InstrumentHeader/ils_relative_response"][:,:,fp,band_idx] # automatically 1/µm
        )

        disp_coef_dict[fp] = l1b_h5["InstrumentHeader/dispersion_coef_samp"][:,fp,band_idx] * u"µm"
        snr_coef_dict[fp] = l1b_h5["InstrumentHeader/snr_coef"][:,:,fp,band_idx]

    end

    @info "Reading retrieval inputs"

    @info "Loading measurement arrays."
    band_name = RE.OCO_band_index_to_string(band_idx)
    meas_full_array = l1b_h5["/SoundingMeasurements/radiance_$(band_name)"][:,:,:]

    if length(size(l1b_h5["/FootprintGeometry/footprint_stokes_coefficients"])) == 4
        stokes_coefs_full_array = l1b_h5["/FootprintGeometry/footprint_stokes_coefficients"][:,:,:,:]
    else
        stokes_coefs_full_array = l1b_h5["/FootprintGeometry/footprint_stokes_coefficients"][1,:,:,:,:]
    end

    sza = l1b_h5["/SoundingGeometry/sounding_solar_zenith"][:,:] * u"°"
    saa = l1b_h5["/SoundingGeometry/sounding_solar_azimuth"][:,:] * u"°"
    relative_velocity = l1b_h5["/SoundingGeometry/sounding_relative_velocity"][:,:] * u"m/s"
    solar_relative_velocity = l1b_h5["/SoundingGeometry/sounding_solar_relative_velocity"][:,:] * u"m/s"
    solar_distance = l1b_h5["/SoundingGeometry/sounding_solar_distance"][:,:] * u"m"
    sounding_qual = l1b_h5["/SoundingGeometry/sounding_qual_flag"][:,:]

    @info "Loading bad samples and spike flags."
    bad_samples = l1b_h5["/InstrumentHeader/bad_sample_list"][:,:, band_idx]
    spike_eof = l1b_h5["/SpikeEOF/spike_eof_weighted_residual_$(band_name)"][:,:,:]

    @info "Done."

    @info "Packing into Dicts."
    for snid in sounding_id_list

        footprint_idx, frame_idx = fp_frame_idx_dict[snid]

        meas_full_dict[snid] = @views meas_full_array[:, footprint_idx, frame_idx]
        stokes_coefs_dict[snid] = @views stokes_coefs_full_array[:, band_idx, footprint_idx, frame_idx]
        sza_dict[snid] = sza[footprint_idx, frame_idx]
        saa_dict[snid] = saa[footprint_idx, frame_idx]
        doppler_factor_dict[snid] = relative_velocity[footprint_idx, frame_idx] ./ RE.SPEED_OF_LIGHT
        solar_doppler_factor_dict[snid] = solar_relative_velocity[footprint_idx, frame_idx] ./ RE.SPEED_OF_LIGHT
        solar_distance_dict[snid] = solar_distance[footprint_idx, frame_idx] ./ RE.AU
        sounding_qual_dict[snid] = sounding_qual[footprint_idx, frame_idx]

    end

    @info "Closing L1B and L2Met files."
    close(l1b_h5)
    close(met_h5)
    if !isnothing(cpr_h5)
        close(cpr_h5)
    end


    return_dict = Dict(
        "fp_frame_idx" => fp_frame_idx_dict,
        "met" => met_dicts,
        "loc" => loc_dict,
        "time" => time_dict,
        "observer" => observer_dict,
        "measured_radiance" => meas_full_dict,
        "stokes_coefs" => stokes_coefs_dict,
        "isrf" => isrf_table_dict,
        "dispersion" => disp_coef_dict,
        "snr_coef" => snr_coef_dict,
        "spike_eof" => spike_eof,
        "bad_samples" => bad_samples,
        "sza" => sza_dict,
        "saa" => saa_dict,
        "doppler_factor" => doppler_factor_dict,
        "solar_doppler_factor" => solar_doppler_factor_dict,
        "solar_distance" => solar_distance_dict,
        "sounding_qual_flag" => sounding_qual_dict
    )

    if !isnothing(cpr_h5)
        return_dict["co2_prior"] = co2_prior_dict
    end

    return return_dict

end


"""
$(TYPEDSIGNATURES)

Reads an aerosol from the NASA ACOS supplementary file, usually named `l2_aerosol_combined.h5`,
and returns an G3RT aerosol property object of type `MieMomAerosolProperty`.

"""
function read_ACOS_aerosol(fname_or_h5::Union{String, Dict}, name::String; wavelength::Bool=true)

    # Open the HDF5 file
    if fname_or_h5 isa String
        h5 = h5open(fname_or_h5, "r")

        # Check if the group even exists
        if !(haskey(h5, name))
            @error "Aerosol $(name) not found."
            close(h5)
            return false
        end

    elseif fname_or_h5 isa Dict
        h5 = fname_or_h5
    end

    # We need to have the root-slash for this to work on dicts as well
    if name[1] != "/"
        name = "/" * name
    end

    if wavelength
        # Flip order if wavelength
        srt = length(h5["$(name)/Properties/wave_number"]):-1:1
        ww = 1e4 ./ h5["$(name)/Properties/wave_number"][:][srt]
        ww_unit = u"µm"
    else
        # Order stays the same
        srt = 1:length(h5["$(name)/Properties/wave_number"])
        ww = h5["$(name)/Properties/wave_number"][:][srt]
        ww_unit = "cm^-1"
    end

    # Conversion to Float64 is needed here since not all types inside this file
    # are always Float64. This is a bit of a hack and might need re-visiting if
    # things break.
    qext = convert(Vector{Float64}, h5["$(name)/Properties/extinction_coefficient"][:][srt])
    qsca = convert(Vector{Float64}, h5["$(name)/Properties/scattering_coefficient"][:][srt])
    ssa = qsca ./ qext

    # Raw ordering in the HDF5 files is: (element, coefficient, wavelength or wavenumber)
    pmom_raw = h5["$(name)/Properties/phase_function_moment"][:,:,:][:,:,srt]
    # Re-shuffle dimensions to be (coefficient, element, wavelength or wavenumber)
    pmom = permutedims(pmom_raw, (2, 1, 3))

    # Reverse sign for γ and ε to match the input convention for XRTM!
    @turbo pmom[:,5,:] *= -1
    @turbo pmom[:,6,:] *= -1


    if h5 isa HDF5.File
        close(h5)
    end

    miemom = RE.MieMomAerosolProperty(
        ww, ww_unit,
        qsca, qext, ssa,
        size(pmom, 1),
        pmom,
        zeros(size(pmom)),
        (:a1, :a2, :a3, :a4, :b1, :b2)
    )

    return miemom

end
