function collect_results(
        s::RE.AbstractSolver,
        buf::RE.EarthAtmosphereBuffer,
)

    # Calculate OE stuff
    q = RE.calculate_OE_quantities(s)
    h = RE.create_pressure_weights(buf.scene.atmosphere)

    # Result dictionary
    result = Dict{String, Any}()

    # Retrieval geometry

    result["/RetrievalGeometry/retrieval_azimuth"] = [buf.scene.observer.viewing_azimuth]
    result["/RetrievalGeometry/retrieval_zenith"] = [buf.scene.observer.viewing_zenith]
    result["/RetrievalGeometry/retrieval_solar_azimuth"] = [buf.scene.solar_azimuth]
    result["/RetrievalGeometry/retrieval_solar_zenith"] = [buf.scene.solar_zenith]


    # State vector results
    sv = s.state_vector
    result["/RetrievedStateVector/state_vector_names"] = permutedims(cat(RE.get_name(sv)..., dims=2), (2,1))
    result["/RetrievedStateVector/state_vector_apriori"] = cat(RE.get_prior_value(sv)..., dims=2)'
    result["/RetrievedStateVector/state_vector_apriori_uncert"] = cat(RE.get_prior_uncertainty(sv)..., dims=2)'
    result["/RetrievedStateVector/state_vector_result"] = cat(RE.get_current_value(sv)..., dims=2)'
    result["/RetrievedStateVector/state_vector_initial"] = cat(RE.get_first_guess(sv)..., dims=2)'
    result["/RetrievedStateVector/state_vector_aposteriori_uncert"] = cat(q.SV_ucert..., dims=2)'

    # Prior covariance matrix
    Sa = s.prior_covariance
    result["/RetrievalResults/apriori_covariance_matrix"] =
        permutedims(Sa[[CartesianIndex()], :, :], (3,2,1))

    # Take the CO2 profile
    gas_co2 = RE.get_gas_from_name(buf.scene.atmosphere, "CO2")
    # If it is present, grab the CO2 profile from the state vector!
    if !isnothing(gas_co2)
        co2_idx = RE.idx_for_profile_sve(gas_co2, s.state_vector)
        co2_profile = RE.get_current_value_with_unit(s.state_vector)[co2_idx] .|> NoUnits
        co2_profile_ap = RE.get_prior_value_with_unit(s.state_vector)[co2_idx] .|> NoUnits
        result["/RetrievalResults/co2_profile"] = cat(co2_profile..., dims=2)'
        result["/RetrievalResults/co2_profile_apriori"] = cat(co2_profile_ap..., dims=2)'
    end

    # Take all XGases
    xgas = RE.calculate_xgas(buf.scene.atmosphere)
    # Store XCO2
    if "CO2" in keys(xgas)
        result["/RetrievalResults/xco2"] = [xgas["CO2"] |> NoUnits]
        # XCO2 uncertainty

    end

    # Measured and modelled, along with measurement noise
    # (as big arrays)
    measured = RE.get_measured(s)
    modeled = RE.get_modeled(s)
    noise = RE.get_noise(s)

    result["/SpectralParameters/modeled_radiance"] = cat(modeled..., dims=2)'
    result["/SpectralParameters/measured_radiance"] = cat(measured..., dims=2)'
    result["/SpectralParameters/measured_radiance_uncert"] = cat(noise..., dims=2)'

    # Calculate residual RMS
    for swin in keys(s.indices)
        meas = RE.get_measured(s, swin, view=false)
        conv = RE.get_modelled(s, swin, view=false)

        rms = sqrt(mean((conv .- meas).^2))
        result["/SpectralParameters/residual_mean_square_$(swin.window_name)"] = [rms]
    end

    # Calculate chi2
    chi2_all = RE.calculate_chi2(s)
    for (swin, chi2) in chi2_all
        result["/SpectralParameters/reduced_chi_squared_$(swin.window_name)"] = [chi2]
    end


    # Store Jacobian
    K = RE.create_K_from_solver(s);
    result["/RetrievalResults/jacobian"] = permutedims(K[[CartesianIndex()], :, :], (3,2,1))

    # High-resolution radiances and Jacobians
    #=
    for (swin, rt) in buf.rt
        result["/HighResSpectra/radiance_$(swin.window_name)"] = rt.hires_radiance.S[:,:]
        for (idx, sve) in enumerate(s.state_vector.state_vector_elements)
            result["/HighResSpectra/jacobian_$(swin.window_name)_$(idx)"] = rt.hires_jacobians[sve].S[:,:]
        end
    end
    =#
    # Normed AK

    # Do this only if we have CO2 in the atmosphere, and are also retrieving it!
    if ("CO2" in keys(xgas)) & !isnothing(q)
        ind = RE.idx_for_profile_sve(gas_co2, s.state_vector)
        if length(ind) > 0
            ak = q.AK[ind,ind];
            ak_norm = (h' * ak)' ./ h
            result["/RetrievalResults/xco2_avg_kernel_norm"] = cat(ak_norm..., dims=2)'
        end
    end

    return result

end
