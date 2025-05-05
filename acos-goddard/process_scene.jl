function process_snid(
    swins,
    scene_inputs,
    state_vector,
    snid,
    buf,
    oe_buf;
    max_iter=10,
    gamma=10.0,
    dsigma_scale=2.0,
    co2_cov=nothing,
    high_options=nothing,
    nus_dict=nothing,
    )


    # Use this for any common quantities
    # (met, etc.)
    common_inputs = first(values(scene_inputs))


    @info "Processing $(snid)"

    N_sv = length(state_vector)
    fp_idx, frame_idx = common_inputs["fp_frame_idx"][snid]

    met = common_inputs["met"]

    # Return if the quality flag signals bad MET data..
    if met[snid]["qflag_met"] != 0
        return nothing
    end


    if haskey(common_inputs, "co2_prior")
        cpr = common_inputs["co2_prior"]
    end

    loc = common_inputs["loc"]
    time = common_inputs["time"]
    observer = common_inputs["observer"]

    doppler_factor = common_inputs["doppler_factor"][snid]
    solar_doppler_factor = common_inputs["solar_doppler_factor"][snid]
    solar_distance = common_inputs["solar_distance"][snid]

    ########################
    # Short-cuts to buffers
    rt_buf = buf.rt_buf
    inst_buf = buf.inst_buf
    ########################

    # These dictionaries map spectrometer index to some quantities, e.g.
    # 1 => stokes coefficient (for spectrometer 1)
    # Per-band Stokes coefficients
    stokes_coefs = Dict(k => v["stokes_coefs"][snid] for (k,v) in scene_inputs)
    # Arrays with measured radiances
    meas = Dict(k => v["measured_radiance"][snid] for (k,v) in scene_inputs)
    # Dispersions
    disps = Dict(k => buf.rt_buf.dispersion[v] for (k,v) in swins)
    # ISRF tables
    isrfs = Dict(k => scene_inputs[k]["isrf"][fp_idx] for (k,v) in swins)

    # Skip any scenes with invalid radiance values
    for (spec, measured) in meas
        if any(measured .< 0)
            @warn "Invalid radiance values found!"
            return nothing
        end
    end

    # Plug in the dispersion prior into the dispersion objects
    for (spec, disp) in disps
        # Empty out first
        empty!(disp.coefficients)
        # Append the dispersion coefficients from the L1b file
        append!(
            disp.coefficients,
            ustrip.(Ref(disp.wavelength_unit), scene_inputs[spec]["dispersion"][fp_idx])
        )
        # Update the object with the new coefficients!
        RE.update_dispersion!(disp)
    end


    # Set the scene location and time
    buf.scene.location = loc[snid]
    buf.scene.time = time[snid]

    # Set the retrieval grid via the ACOS method
    @views buf.scene.atmosphere.pressure_levels[:] = ustrip(
            RE.create_ACOS_pressure_grid(met[snid]["surface_pressure_met"]) .|>
            buf.scene.atmosphere.pressure_unit
        )

    #=
        The meteorological quantities (p,T,q)
    =#

    # Pressure grid must be turned into Pa
    RE.ingest!(
            buf.scene.atmosphere,
            :met_pressure_levels,
            met[snid]["pressure_levels_met"]
            )


    # Ingest specific humidity, which has unit 1 (kg/kg)
    RE.ingest!(
        buf.scene.atmosphere,
        :specific_humidity_levels,
        met[snid]["specific_humidity_profile_met"]
    )

    # Temperature profile must be turned into K
    RE.ingest!(
        buf.scene.atmosphere,
        :temperature_levels,
        met[snid]["temperature_profile_met"]
    )

    #=
        Supply the CO2 prior profile
    =#
    # Grab the CO2 atmosphere element, if it exists
    atm_co2 = RE.get_gas_from_name(buf.scene.atmosphere, "CO2")

    if !isnothing(atm_co2)

        co2_met_levels = ustrip.(atm_co2.vmr_unit, cpr[snid]["co2_prior"])

        # Interpolate onto the retrieval grid (this needs to be better! it's not Xgas conserving!)
        co2_prior = RE.atmospheric_profile_interpolator_linear(
            buf.scene.atmosphere.met_pressure_levels,
            co2_met_levels,
            buf.scene.atmosphere.pressure_levels
            )

        # Move the CO2 profile into the CO2 gas object
        @views atm_co2.vmr_levels[:] .= co2_prior
    end

    # From the MET q (spec. humidity) profile, we interpolate down to the retrieval grid,
    # so we can use it as the H2O VMR profile. Of course we only do that if a water vapor
    # gas object is present..

    atm_h2o = RE.get_gas_from_name(buf.scene.atmosphere, "H2O")

    if !isnothing(atm_h2o)

        sh_profile = RE.atmospheric_profile_interpolator_linear(
            buf.scene.atmosphere.met_pressure_levels,
            buf.scene.atmosphere.specific_humidity_levels,
            buf.scene.atmosphere.pressure_levels
            )

        h2o_profile = @. sh_profile / (1.0 - sh_profile) * RE.MM_AIR_TO_H2O
        @views atm_h2o.vmr_levels[:] .= ustrip.(Ref(atm_h2o.vmr_unit), h2o_profile*Unitful.NoUnits)

    end

    # Set oxygen VMR if an O2 gas object is present
    atm_o2 = RE.get_gas_from_name(buf.scene.atmosphere, "O2")
    if !isnothing(atm_o2)
        atm_o2.vmr_levels[:] .= 0.2095
    end

    # Calculate altitude and gravity
    RE.calculate_altitude_and_gravity!(buf.scene)
    # Calculate mid-layer values
    RE.calculate_layers!(buf.scene.atmosphere)

    # Ingest the observer object into the scene, along with solar angles


    # IMPORTANT!!!!
    # (taken from L2FP)
    #=
      Azimuth is modified because the convention used by the OCO L1B file is to
      take both solar and observation angles as viewed from an observer
      standing in the FOV.  In this convention, the convention for glint
      would be a relative azimuth difference of 180 degrees, as the
      spacecraft and sun would be on opposite sides of the sky. However, the
      radiative transfer convention is that the azimuth angles must be the
      same for glint (it is "follow the photons" convention). However, we'd
      like the solar azimuth to not be changed, so as to continue to agree
      with zenith, so this change of the observation azimuth has the effect
      of putting everything in a "reverse follow-the-photons" convention,
      where we look from the satellite to the FOV, then from the FOV to the
      sun.  Note that because of an old historical reason, however, both
      zenith angles remain > 0 and < 90, even in the RT convention.
    =#

    buf.scene.observer = observer[snid]

    #saa = common_inputs["saa"][snid] |> u"°" |> ustrip
    #vaa = buf.scene.observer.viewing_azimuth
    #buf.scene.observer.viewing_azimuth = mod(vaa + 180 - saa, 360)

    @info "Observer viewing zenith: $(buf.scene.observer.viewing_zenith)"
    @info "Observer viewing azimuth: $(buf.scene.observer.viewing_azimuth)"
    buf.scene.solar_zenith = common_inputs["sza"][snid] |> u"°" |> ustrip # in degrees
    buf.scene.solar_azimuth = common_inputs["saa"][snid] |> u"°" |> ustrip # in degrees
    @info "Solar zenith: $(buf.scene.solar_zenith)"
    @info "Solar azimuth: $(buf.scene.solar_azimuth)"

    #=
    ##################################################
    # Prepare the SV with priors and prior covariances
    ##################################################
    =#


    # Aerosols!

    # AOD retrieval
    for (idx, sve) in RE.StateVectorIterator(state_vector, RE.AerosolOpticalDepthSVE)
        if sve.aerosol.aerosol_name == "ST" # Stratospheric aerosols
            sve.prior_value = -5.116
            sve.first_guess = -5.116
            sve.prior_covariance = (1.8)^2
        elseif sve.aerosol.aerosol_name == "Water" # Water cloud
            sve.prior_value = -4.382
            sve.first_guess = -4.382
            sve.prior_covariance = (1.8)^2
        elseif sve.aerosol.aerosol_name == "Ice" # Ice cloud
            sve.prior_value = -4.382
            sve.first_guess = -4.382
            sve.prior_covariance = (1.8)^2
        else # Tropospheric type!
            sve.prior_value = log(sve.aerosol.total_optical_depth)
            sve.first_guess = log(sve.aerosol.total_optical_depth)
            sve.prior_covariance = (0.5)^2
        end
    end

    # height retrieval
    for (idx, sve) in RE.StateVectorIterator(state_vector, RE.AerosolHeightSVE)
        if sve.aerosol.aerosol_name == "ST" # Stratospheric aerosols
            sve.prior_value = 0.03
            sve.first_guess = 0.03
            sve.prior_covariance = (1.0e-4)^2
        elseif sve.aerosol.aerosol_name == "Water" # Water cloud
            sve.prior_value = 0.75
            sve.first_guess = 0.75
            sve.prior_covariance = (0.4)^2
        elseif sve.aerosol.aerosol_name == "Ice" # Ice cloud
            # Ice cloud height is 100 hPa below the tropopause height
            ptropo = met[snid]["tropopause_pressure_met"]
            ice_height = ustrip(buf.scene.atmosphere.pressure_unit, ptropo - 100.0u"hPa") /
                buf.scene.atmosphere.pressure_levels[end]
            @info "Tropopause height is: $(ptropo)"
            @info "Setting ice cloud height to $(ice_height)."
            sve.prior_value = ice_height
            sve.first_guess = ice_height
            sve.prior_covariance = (0.2)^2
        else # Tropospheric type!
            sve.prior_value = 0.9
            sve.first_guess = 0.9
            sve.prior_covariance = (0.2)^2
        end
    end

    # width retrieval
    for (idx, sve) in RE.StateVectorIterator(state_vector, RE.AerosolWidthSVE)
        if sve.aerosol.aerosol_name == "ST" # Stratospheric aerosols
            sve.prior_value = 0.04
            sve.first_guess = 0.04
            sve.prior_covariance = (0.01)^2
        elseif sve.aerosol.aerosol_name == "Water" # Water cloud
            sve.prior_value = 0.1
            sve.first_guess = 0.1
            sve.prior_covariance = (0.01)^2
        elseif sve.aerosol.aerosol_name == "Ice" # Ice cloud
            sve.prior_value = 0.04
            sve.first_guess = 0.04
            sve.prior_covariance = (0.01)^2
        else # Tropospheric type!
            sve.prior_value = 0.05
            sve.first_guess = 0.05
            sve.prior_covariance = (0.01)^2
        end
    end

    # Ingest dispersion polynomial coefficient priors
    for (idx, sve) in RE.StateVectorIterator(state_vector, RE.DispersionPolynomialSVE)

        # What is the coefficient order of this SVE?
        o = sve.coefficient_order
        # To which spectral window is this SVE attached to?
        disp_swin = sve.dispersion.spectral_window
        # What is the spectrometer index that belongs to this swin?
        # (note that this might break in pathological situations like double spectral windows..)
        spec = first(k for (k,v) in swins if v == disp_swin)

        sve.first_guess = ustrip(sve.unit, scene_inputs[spec]["dispersion"][fp_idx][o+1])
        sve.prior_value = ustrip(sve.unit, scene_inputs[spec]["dispersion"][fp_idx][o+1])
        # These need to be rather tight, the Jacobian for dispersion polynomial coefficients
        # is quasi-valid for a very small range only. If the initial dispersion is too far
        # away from a good value, the retrieval will likely struggle to produce anything
        # meaningful anyway.
        sve.prior_covariance = (1e-5 * sve.prior_value)^2

    end

    # Set CO2 state vector
    for (idx, sve) in RE.StateVectorIterator(state_vector, RE.GasVMRProfileSVE)

        # Make sure it's our desired one
        if sve.gas === RE.get_gas_from_name(buf.scene.atmosphere, "CO2")
            sve.first_guess = sve.gas.vmr_levels[sve.level]
            sve.prior_value = sve.first_guess
            # Prior covariance will be set later on..
            sve.prior_covariance = ustrip(sve.unit, 0.0u"ppm")
        end

    end


    # Set psurf
    for (idx, sve) in RE.StateVectorIterator(state_vector, RE.SurfacePressureSVE)

        sve.prior_value = ustrip(sve.unit,
            buf.scene.atmosphere.pressure_levels[end] * buf.scene.atmosphere.pressure_unit)
        sve.first_guess = sve.prior_value

    end

    # Set a Band 1 ZLO, and restrict to ~1% of the continuum
    for (idx, sve) in RE.StateVectorIterator(state_vector, RE.ZeroLevelOffsetPolynomialSVE)
        if (sve.swin === swins[1]) & (sve.coefficient_order == 0)
            # Make this only work for spectral band 1
            cont_level = percentile(meas[1], 95)
            zlo_prior_std = 0.01 * cont_level
            @info "Setting band 1 ZLO prior covariance to $(zlo_prior_std)."
            sve.prior_covariance = (zlo_prior_std)^2
        end
    end


    # Calculate prior albedo value from radiances, and also
    # account for the Stokes coefficient
    for (spec, swin) in swins
        # Find out any unit conversion factor between measured radiance (rt_buf.radiance_unit)
        # and the radiance we are using internally.
        unit_fac = rt_buf.radiance_unit / buf.rt[swin].radiance_unit

        # Calculate apparent albedo from the measured radiances

        signal = percentile(meas[spec], 99)
        if buf.rt[swin].solar_model isa RE.OCOHDFSolarModel
            solar_strength = maximum(buf.rt[swin].solar_model.continuum) * stokes_coefs[spec][1]
        else
            solar_idx = searchsortedfirst.(Ref(buf.rt[swin].solar_model.ww), swin.ww_grid)
            solar_strength = maximum(buf.rt[swin].solar_model.irradiance[solar_idx]) * stokes_coefs[spec][1]
        end

        albedo_prior = pi * signal / (solar_strength * cosd(buf.scene.solar_zenith)) * unit_fac

        # Ingest prior albedo into state vector values
        for (idx, sve) in RE.StateVectorIterator(state_vector, RE.BRDFPolynomialSVE)
            if sve.swin === swin
                if sve.coefficient_order == 0
                    sve.prior_value = albedo_prior
                    sve.first_guess = albedo_prior
                    sve.prior_covariance = 1.0
                else
                    sve.prior_value = 0
                    sve.first_guess = 0
                    sve.prior_covariance = 10.0
                end
            end
        end
    end

    # Empty out all state vector element iterations,
    # and populate with first guesses
    for sve in state_vector.state_vector_elements
        empty!(sve.iterations)
        push!(sve.iterations, sve.first_guess)
    end

    # Define the boundaries for LSI here:
    LSI_bounds = Dict{Int, Vector{Float64}}()
    if 1 in keys(swins)
        LSI_bounds[1] = [
            0, 0.001, 0.01, 0.05, 0.1, 0.15, 0.225, 0.338, 0.508,
            0.7, 0.85, 1, 1.15, 1.4, 1.719, 2.581, 3.875, 5.817,
            8.733, 13, 20, 30, 44, 70, 100, 200, 300, 1e+06
        ]
    end
    if 2 in keys(swins)
        LSI_bounds[2] = [
            0, 0.00027, 0.001, 0.004, 0.011, 0.025, 0.05, 0.1,
            0.2, 0.4, 0.7, 1, 1.25, 1.5, 2, 1e+06
            ]
    end
    if 3 in keys(swins)
        LSI_bounds[3] = [
            0, 0.011, 0.03, 0.1, 0.3, 0.5, 0.85, 1.03, 1.21, 1.4, 1.8, 2.3, 3,
            3.7, 4.5, 5.7, 7, 9.5, 12.5, 17, 25, 1e+06
            ]
    end


    # Forward model keyword arguments, which are needed to
    # be passed to the solver.

    fm_kwargs = (
        buf=buf,
        inst_buf=inst_buf,
        oe_buf=oe_buf,
        rt_buf=rt_buf,
        dispersion=disps,
        isrf=isrfs,
        stokes=stokes_coefs,
        NUS=nus_dict,
        LSI_bounds=LSI_bounds,
        doppler_factor=doppler_factor,
        solar_doppler_factor=solar_doppler_factor,
        solar_distance=solar_distance,
        high_options=high_options
    )

    # Create a view to the Sa buffer, such that the
    # solver can use it later on.
    Sa = @views oe_buf.Sa[1:N_sv, 1:N_sv]

    # Calculate measurement noise
    # NOTE:
    # Some of the inverse method linear algebra might
    # produce NaNs due to Float32 issues, so I am forcing
    # this particular array to be a Float64.

    MaxMS = [7.0e20, 2.45e20, 1.25e20]
    meas_noise = Dict{Int, Vector{Float64}}()
    for spec in keys(swins)
        # Calculate the noise-equivalent radiance according to ACOS
        meas_noise[spec] = convert.(Float64, RE.OCO_calculate_noise(
                meas[spec],
                common_inputs["snr_coef"][fp_idx],
                MaxMS[spec]
                )
            )

        # Check the spike flag and enhance noise when spikes occur
        loc_spike = findall(scene_inputs[spec]["spike_eof"][:,fp_idx,frame_idx] .>= 3)
        meas_noise[spec][loc_spike] .= MaxMS[spec]

    end

    # Create a solver object for our chosen solver
    solver = RE.LMSolver(
        forward_model!,
        state_vector,
        Sa,
        gamma,
        max_iter,
        dsigma_scale,
        Dict(swins[spec] => disps[spec] for spec in keys(swins)),
        rt_buf.indices,
        rt_buf.radiance,
        rt_buf.jacobians,
        Dict(disps[spec] => meas[spec] for spec in keys(swins)),
        Dict(disps[spec] => meas_noise[spec] for spec in keys(swins))
        )


    # Adjust the prior covariance if wanted
    if !isnothing(co2_cov)
        # CO2 cov is a matrix here
        for (idx1,sve1) in RE.StateVectorIterator(state_vector, RE.GasVMRProfileSVE)
            for (idx2,sve2) in RE.StateVectorIterator(state_vector, RE.GasVMRProfileSVE)
                solver.prior_covariance[idx1, idx2] = co2_cov[sve1.level, sve2.level]
            end
        end
    end

    # Perform the iterations
    iter_count = 0
    converged = false

    while (iter_count < solver.max_iterations)

        converged = RE.check_convergence(solver, verbose=true)
        if converged
            @debug "Reached convergence!"
            break
        end

        @time iter_success = RE.next_iteration!(solver; fm_kwargs)
        iter_count += 1

        if !iter_success
            @debug "FAILED iteration #$(iter_count)"
            break
        else
            @info "Successful iteration #$(iter_count)"
            @info "Length of SV: $(RE.get_iteration_count(solver))"
        end

        chi2 = RE.calculate_chi2(solver)
        for (k,v) in chi2
            @info "$(k): χ² =  $(v)"
        end

    end
    @debug "Iterations done."

    # Update the atmosphere one final time!
    # This pushes updates of the SV back into the atmosphere
    # (necessary for other routines that calculate
    #  quantities from the atmosphere rather than the SV)

    # Update atmosphere
    RE.atmosphere_statevector_update!(
        buf.scene.atmosphere,
        state_vector
    )

    # Update atmosphere elements
    RE.atmosphere_element_statevector_update!(
        buf.scene.atmosphere.atm_elements,
        state_vector
    )

    q = RE.calculate_OE_quantities(solver);
    if !isnothing(q)
        RE.print_posterior(solver, q)
    end

    return solver, fm_kwargs

end
