"""

## Details
`nus` is some vector (with the same length as `swin.ww_grid` that evaluates to `true` for all
spectral points for which valid hires RT calculations have been performed. Conversely,
`nus` must evaluate to false for those points for which this function will provide interpolated
values.

"""
function NUS_correction!(
        rt::RE.MonochromaticRTMethod,
        swin::RE.AbstractSpectralWindow,
        nus::AbstractVector;
        kind="linear"
    )

    knots_idx = findall(nus)
    fills_idx = findall(.!nus)
    N_fills = length(fills_idx)

    # We interpolate in wavelength or wavenumber space
    knots = swin.ww_grid[knots_idx]
    fills = swin.ww_grid[fills_idx]

    # Loop over Stokes components `s`
    for s in 1:size(rt.hires_radiance, 2)

        # Do a quick check: all radiance values that are in `fills` should
        # have some non-zero value.
        N_zero = (rt.hires_radiance[knots_idx, s] .== 0) |> sum
        N_zero > 0 && @warn "NUS problem: N=$(N_zero) zeros found in radiance (stokes=$(s))."

        # Build interpolation object for radiance
        itp = linear_interpolation(
            knots, rt.hires_radiance.S[knots_idx, s],
            extrapolation_bc=Flat()
            )

        # Interpolate
        for i in 1:N_fills
            rt.hires_radiance.S[fills_idx[i], s] = itp(fills[i])
        end

        # Do the same for each used weighting function
        for (k, wf_idx_list) in rt.wfunctions_map
            for wf in wf_idx_list

                this_wf = @view rt.hires_wfunctions[wf].S[:, s][knots_idx]

                itp_wf = linear_interpolation(
                    knots, this_wf,
                    extrapolation_bc=Flat()
                    )

                # Interpolate at `fills` locations (or `fills_idx` for the array positions to be filled)
                for i in 1:N_fills
                    rt.hires_wfunctions[wf].S[fills_idx[i], s] = itp_wf(fills[i])
                end
            end
        end

    end


    #=
        And finally re-calculate all Jacobians using the reconstructed weighting functions.

        THIS IS IMPORTANT. Though we reconstructed the full weighting functions, the Jacobians
        are *not* corrected at this point, so they would contain a lot of zeros for the spectral
        points which are not evaluated during the monochromatic RT loop.

        This might cost a little bit of computational effort, but for e.g. OCO-type retrievals,
        the overall time for this is quite negligible (fractions of a second).

    =#
    for (i, sve) in enumerate(rt.state_vector.state_vector_elements)
        if RE.calculate_jacobian_before_isrf(sve)
            RE.calculate_rt_jacobian!(rt.hires_jacobians[sve], rt, sve)
        end
    end


end