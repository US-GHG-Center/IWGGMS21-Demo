"""
Creates a statevector for ACOS OCO retrievals.

# Details
Note that many of the state vector values (priors, first guess, prior covariance etc., are
*placeholders* only, as they will be changed later on when the scene-dependent configuration
is performed. We need the state vector, however, to create the buffer object.
"""
function create_statevector(
        swins::Dict,
        disps::Dict,
        atm_elements::Vector{<:RE.AbstractAtmosphereElement}
    )


    # Get the array of present spectrometers
    specs = collect(keys(swins))
    # Might as well sort..
    sort!(specs)

    # Allocate empty vector of state vector elements
    sv = RE.AbstractStateVectorElement[]

    @info "Creating the state vector!"
    @info "..."

    # .. now populate depending on the choice of spectrometers

    # Retrieve surface pressure only if O2 A-band is present!
    if 1 in specs
        @info "Retrieving surface pressure."

        sve_psurf = RE.SurfacePressureSVE(
            u"hPa", # Choose a unit here, we can make this any pressure unit
            0.0, # This will be filled in later
            0.0, # So will this
            (4.0)^2, # Choose an appropriate prior covariance
            )

        push!(sv, sve_psurf)
    end


    # Gaussian aerosol AOD and height retrieval

    for atm in atm_elements
        if atm isa RE.GaussAerosol

            @info "Retrieving log-AOD for $(atm)."
            sve_aod = RE.AerosolOpticalDepthSVE(
                atm, # reference to the aerosol object itself
                true, # Log-scale?
                Unitful.NoUnits,
                0.0, # prior
                0.0, # first guess
                0.0, # prior covariance
                )

            push!(sv, sve_aod)

            @info "Retrieving fractional height for $(atm)."
            sve_height = RE.AerosolHeightSVE(
                atm, # reference to the aerosol object itself
                false, # Log-scale?
                Unitful.NoUnits,
                0.0, # prior
                0.0, # first guess
                0.0, # prior covariance
                )

            push!(sv, sve_height)

            @info "Retrieving fractional width for $(atm)"
            sve_width = RE.AerosolWidthSVE(
                atm, # reference to the aerosol object itself
                false, # Log-scale?
                Unitful.NoUnits,
                0.0, # prior
                0.0, # first guess
                0.0, # prior covariance
            )

            push!(sv, sve_width)

        end
    end


    # T offset will stay the same, this one is easy
    @info "Retrieving T profile offset."
    sve_T_offset = RE.TemperatureOffsetSVE(
        u"K",
        0.0, # First guess
        0.0, # Prior value
        (5.0)^2, # Prior covariance
    );
    push!(sv, sve_T_offset)


    # Gas scale retrieval for H2O if spectrometers 2 or 3 are involved
    if (2 in specs) | (3 in specs)

        @info "Retrieving H2O VMR scale factor."
        # Grab the H2O gas object
        gas_h2o = filter(x -> x isa RE.GasAbsorber && x.gas_name == "H2O", atm_elements)[1]

        sve_h2o_scale = RE.GasLevelScalingFactorSVE(
            1, # Start level index
            length(gas_h2o.vmr_levels), # End level index
            gas_h2o, # Gas absorber which this SV element references
            Unitful.NoUnits, # Unit of scaling factor (could be e.g. percent)
            1., # First guess
            1., # Prior value
            (0.5)^2, # Prior covariance
        );

        push!(sv, sve_h2o_scale)

    end


    # CO2 VMR profile retrieval if spectrometers 2 or 3 are involved
    if (2 in specs) | (3 in specs)

        @info "Retrieving CO2 VMR profile."
        # Grab the H2O gas object
        gas_co2 = filter(x -> x isa RE.GasAbsorber && x.gas_name == "CO2", atm_elements)[1]

        # This generates a new set of SVEs, one for each layer, which is why
        # we need the splatting operator `...` to add all of them separately, rather than
        # adding a vector of SVEs.
        sve_co2_profile = RE.GasVMRProfileSVE(
            length(gas_co2.vmr_levels),
            gas_co2, # Gas absorber
            Unitful.ppm, # Unit
        )

        push!(sv, sve_co2_profile...)

        sve_co2_scale = RE.GasLevelScalingFactorSVE(
            1, # Start level index
            20, # End level index
            gas_co2, # Gas absorber which this SV element references
            Unitful.NoUnits, # Unit of scaling factor (could be e.g. percent)
            1.0, # First guess
            1.0, # Prior value,
            (1.0)^2, # Prior covariance
        )

        #push!(sv, sve_co2_scale)

    end


    # Add spectrometer-specific SVEs (surface, dispersion, ZLO, ...)
    for spec_idx in specs

        # Add ZLO to mimic SIF (not yet implemented)
        if spec_idx == 1
            @info "Retrieving ZLO order 0 for band $(spec_idx)"

            sve_zlo_band1 = RE.ZeroLevelOffsetPolynomialSVE(
                swins[spec_idx], # Retrieval window name
                0, # ZLO polynomial order
                u"µm", # wavelength units
                u"ph/s/m^2/sr/µm", #radiance units work
                0.0, # First guess
                0.0, # Prior value
                0.0 # Prior cov
            )

            push!(sv, sve_zlo_band1)

        end

        # Add up to second order BRDF polynomials
        for o in [0,1,2]

            #if (spec_idx == 1) & (o > 1)
                # Second order polynomial and higher really
                # does not work well with the O2 A-band..
            #    continue
            #end

            @info "Retrieving BRDF polynomial order $(o) for band $(spec_idx)."
            sve_brdf = RE.BRDFPolynomialSVE(
                swins[spec_idx],
                RE.LambertianPolynomialKernel,
                o, # Polynomial order
                u"µm", # Associated unit
                0.0, # First guess
                0.0, # Prior value
                0.0, # Prior covariance
            )

            push!(sv, sve_brdf)

        end

        # Add up to linear dispersion polynomials
        for o in [0,1]

            @info "Retrieving dispersion polynomial order $(o) for band $(spec_idx)."
            sve_dispersion = RE.DispersionPolynomialSVE(
                disps[spec_idx],
                o, # Polynomial order
                u"μm", # Unit
                0.0, # First guess
                0.0, # Prior value
                0.0, # Prior covariance
            )

            push!(sv, sve_dispersion)

        end

    end

    # All done, create and return
    return RE.RetrievalStateVector(sv)

end
