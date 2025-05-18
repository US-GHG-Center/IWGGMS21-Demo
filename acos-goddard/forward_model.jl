function forward_model!(
    SV::RE.AbstractStateVector;
    buf::RE.EarthAtmosphereBuffer,
    inst_buf::RE.InstrumentBuffer,
    oe_buf::RE.OEBuffer,
    rt_buf::RE.AbstractRTBuffer,
    dispersion::Dict,
    isrf::Dict,
    stokes::Dict,
    NUS::Union{Nothing, Dict},
    LSI_bounds::Union{Nothing, Dict},
    doppler_factor::AbstractFloat,
    solar_doppler_factor::AbstractFloat,
    solar_distance::AbstractFloat, # in multiples of AU
    high_options::Union{Nothing, <:AbstractDict, Vector}
    )

    spectral_windows = Dict(spec => disp.spectral_window for (spec, disp) in dispersion)

    # Update the dispersion object based on the state vector
    # (You must make sure that the dispersion state vector elements
    #  are bound to the used dispersion!)
    for (spec, disp) in dispersion
        RE.update_dispersion!(disp, SV)
    end

    # Re-calculate indices for all windows
    RE.calculate_indices!(buf)

    # Update the solar scaler field in the RT container
    for (swin, rt) in buf.rt
        RE.solar_scaler_statevector_update!(rt, SV)
    end

    # If we have a surface SV, plug in the current value
    for (swin, rt) in buf.rt
        RE.surfaces_statevector_update!(rt.scene, SV)
    end

    #RE.update_specific_humidity_from_H2O!(buf.scene.atmosphere)

    # If we retrieve surface pressure, create a new grid here!
    # For ACOS, we don't just shift the surface pressure level, but create
    # a new grid altogether, based on the new surface pressure.

    for (idx, sve) in RE.StateVectorIterator(SV, RE.SurfacePressureSVE)

        # Create new pressure grid!
        new_psurf = RE.get_current_value(sve)
        @views buf.scene.atmosphere.pressure_levels[:] = ustrip.(
            Ref(buf.scene.atmosphere.pressure_unit),
            RE.create_ACOS_pressure_grid(new_psurf * sve.unit)
            )


        # Shift lowest pressure level according to the psurf SVE
        # (and account for units)
        #=
        new_psurf = ustrip(
            buf.scene.atmosphere.pressure_unit,
            RE.get_current_value_with_unit(sve)
            )

        buf.scene.atmosphere.pressure_levels[end] = new_psurf
        =#
        # .. and calculate the mid-layer values that need to be
        # inserted into the atmosphere object.
        @views buf.scene.atmosphere.pressure_layers[:] =
            RE.levels_to_layers(buf.scene.atmosphere.pressure_levels)

        # Re-calculate z and g
        RE.calculate_altitude_and_gravity!(buf.scene)

        # .. and again calculate mid-layer values.
        buf.scene.atmosphere.gravity_layers[:] =
            RE.levels_to_layers(buf.scene.atmosphere.gravity_levels)
        buf.scene.atmosphere.altitude_layers[:] =
            RE.levels_to_layers(buf.scene.atmosphere.altitude_levels)

    end

    # Obtain down-sampled solar radiance at the required wavelength grid
    for (swin, rt) in buf.rt

        RE.calculate_solar_irradiance!(
            rt,
            swin,
            rt.solar_model,
            doppler_factor=solar_doppler_factor
        )
        # And scale by this relative solar distance factor
        @views rt.hires_solar.I[:] /= (solar_distance^2)
    end

    # Given the current statevector - modify any atmospheric element
    # according to whatever the state vector element dictates.
    # (e.g. scale the VMR of a gas according to the SV)
    RE.atmosphere_element_statevector_update!(buf.scene.atmosphere.atm_elements, SV)

    # Update the atmosphere according to the current state vector
    RE.atmosphere_statevector_update!(buf.scene.atmosphere, SV)

    # Calculate optical properties inside buffer
    RE.calculate_earth_optical_properties!(buf, SV, N_sublayer=10)

    # Empty out all RT results, since we might want to use additive methods
    # later on.
    for rt in values(buf.rt)
        RE.clear!(rt)
    end

    #=
        Perform RT
    =#

    for (spec, swin) in spectral_windows

        # If non-uniform sampling (NUS) is wanted and present, we must add the NUS array
        # into the "sampling" key of the model option dictionary, so that the RT routines
        # know to skip these points.
        if !isnothing(NUS)
            for opts in buf.rt[swin].model_options
                opts["sampling"] = NUS[spec]
            end
        end

        # Minor hack: if "calc_derivs" is not present in
        # buf.rt[swin].model_options["options"], that's probably a good sign that the user
        # does not want derivatives for LSI either. We remove that option here.
        for lmo in buf.rt[swin].model_options
            if !isnothing(high_options) & !("calc_derivs" in lmo["options"])

                @debug "Removed -calc_derivs- from LSI calculations"
                for mo in high_options
                    filter!(!=("calc_derivs"), mo["options"])
                end
                break # No need to go on
            end
        end



        #####################################################
        # Calculate !!
        #####################################################
        @info "Performing RT calculations for $(swin) ..."
        RE.calculate_radiances_and_jacobians!(buf.rt[swin])

        # Interpolate NUS points
        if !isnothing(NUS)
            NUS_correction!(
                buf.rt[swin],
                swin,
                NUS[spec];
                kind="linear"
                )
        end


        # Peform LSI correction if `high_options` are supplied.
       if !isnothing(high_options)
            @info "(LSI correction ... )"
            # Create the method
            lsi = RE.LSIRTMethod(
                LSI_bounds[spec],
                buf.rt[swin],
                high_options
            )

            # Run the binned RT calculations, and perform the correction itself.
            RE.perform_LSI_correction!(lsi)

        end

    end

    @info " ... done!"

    # Apply Stokes coefficients and solar model to produce the at-instrument
    # radiances *before* applying the ISRF. Note that we store the intermediate
    # radiance result in the rt.optical_properties.Nhi1 temporary vector.
    # Finally, we apply the ISRF to get the proper at-detector result.
    for (spec, swin) in spectral_windows

        # Grab the RT object
        rt = buf.rt[swin]
        Nstokes = size(rt.hires_radiance, 2)

        # Temp object for holding results
        # IMPORTANT! THIS WILL KEEP THE HIRES RADIANCES FOR LATER TOO!
        hires = rt.optical_properties.tmp_Nhi1
        hires[:] .= 0


        @debug "Applying Stokes coeffs for $(swin): $(stokes[spec])"
        # Radiance
        for s in 1:Nstokes
            @turbo hires[:] .+= rt.hires_radiance.S[:,s] .* stokes[spec][s]
        end

        # Multiply by solar model
        @turbo hires[:] .*= rt.hires_solar.S[:,1]


        success = RE.apply_isrf_to_spectrum!(
                inst_buf, # Convolution buffer
                isrf[spec], # ISRF table
                dispersion[spec], # Dispersion
                hires, # radiance
                swin,
                doppler_factor=doppler_factor
            )

        if !success
            @warn "Application of ISRF on radiances for $(swin) FAILED."
            return false
        else
            @debug "Application of ISRF on radiances for $(swin)successful."
        end

        # Store in RT buffer, but account for unit differences!
        @views rt_buf.radiance.I[rt_buf.indices[swin]] =
            inst_buf.low_res_output[dispersion[spec].index] *
            buf.rt[swin].radiance_unit / rt_buf.radiance_unit

    end

    # Apply any radiance correction post-ISRF, such as ZLO
    for sve in SV.state_vector_elements
        RE.apply_radiance_correction!(rt_buf, sve)
    end


    # Same as above, but we do one Jacobian at a time
    for (spec, swin) in spectral_windows

        # Grab the RT object
        rt = buf.rt[swin]
        Nstokes = size(rt.hires_radiance, 2)

        # Temp object for holding results
        hires = rt.optical_properties.tmp_Nhi2

        for (i, sve) in enumerate(SV.state_vector_elements)
            if RE.calculate_jacobian_before_isrf(sve)

                # Zero out
                hires[:] .= 0
                # Apply Stokes coefficient and add to result
                for s in 1:Nstokes
                    @turbo hires[:] .+= rt.hires_jacobians[sve].S[:,s] .* stokes[spec][s]
                end

                # Multiply by solar model
                @turbo hires[:] .*= rt.hires_solar.S[:,1]

                # Apply ISRF
                success = RE.apply_isrf_to_spectrum!(
                        inst_buf, # Convolution buffer
                        isrf[spec], # ISRF table
                        dispersion[spec], # Dispersion
                        hires,
                        swin,
                        doppler_factor=doppler_factor
                    )

                if !success
                    @warn "Application of ISRF on Jacobian $(sve) for $(swin) FAILED."
                    return false
                else
                    @debug "Application of ISRF on Jacobian $(sve) for $(swin) successful."
                end

                # Move over to RT Buffer object
                @views rt_buf.jacobians[sve].I[rt_buf.indices[swin]] =
                    inst_buf.low_res_output[dispersion[spec].index]

            end # End if Jacobian should be calculated BEFORE ISRF
        end # End SV loop
    end # End spectral window loop


    # Post-ISRF Jacobian calculations for ZLO

    for (i, sve) in RE.StateVectorIterator(SV, RE.ZeroLevelOffsetPolynomialSVE)
        RE.calculate_jacobian!(rt_buf, sve)
    end

    # Post-ISRF Jacobian calculation for ISRF
    for (i, sve) in RE.StateVectorIterator(SV, RE.DispersionPolynomialSVE)

        for (spec, swin) in spectral_windows

            if !(sve.dispersion.spectral_window === swin)
                # This is not our spectral window - SKIP!
                continue
            end

            # This has the hires radiances..
            hires_rad = buf.rt[swin].optical_properties.tmp_Nhi1
            success = RE.calculate_dispersion_polynomial_jacobian!(
                inst_buf.low_res_output,
                inst_buf,
                sve,
                isrf[spec],
                dispersion[spec],
                hires_rad,
                swin
            )

            # Re-set to zero! IMPORTANT!
            @views rt_buf.jacobians[sve].S[:] .= 0
            # Ingest jacobian into the right place
            rt_buf.jacobians[sve].I[rt_buf.indices[swin]] =
                inst_buf.low_res_output[dispersion[spec].index]

            if !success
                @warn "Application of ISRF on Jacobian $(sve) for $(swin) FAILED."
                return false
            else
               @debug "Application of ISRF on Jacobian $(sve) for $(swin) successful."
            end

        end

    end

    # Correct *all* Jacobians for unit differences
    for (sve, jac) in rt_buf.jacobians
        for (spec, swin) in spectral_windows
            @views jac.I[rt_buf.indices[swin]] *= 1.0 * buf.rt[swin].radiance_unit / rt_buf.radiance_unit
        end
    end

    # Important!
    # Re-set the atmosphere to its prior state!
    RE.atmosphere_element_statevector_rollback!(buf.scene.atmosphere.atm_elements, SV)
    RE.atmosphere_statevector_rollback!(buf.scene.atmosphere, SV)


    return true

end
