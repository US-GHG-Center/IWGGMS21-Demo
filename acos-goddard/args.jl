function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table! s begin

        "--logfile"
            help = "Path to log file"
            arg_type = String
            required = false
        "--solar_model"
            help = "Path to solar model file"
            arg_type = String
            required = true
        "--L1b"
            help = "Path to L1b file"
            arg_type = String
            required = true
        "--L2Met"
            help = "Path to L2Met file"
            arg_type = String
            required = true
        "--L2CPr"
            help = "Path to L2CPr file"
            arg_type = String
            required = true
        "--sounding_id"
            help = "Sounding ID (16-digit integer)"
            arg_type = Int
            required = false
        "--sounding_id_list"
            help = "Path to Sounding ID list"
            arg_type = String
            required = false
        "--spec"
            help = "Which spectrometers? Separated by comma: e.g. 1,2 or 1, or 1,2,3 (O2, weak CO2, strong CO2)"
            arg_type = String
            required = true
        "--aerosols"
            help = "Whether to use aerosols or not (default = true)"
            arg_type = Bool
            required = false
            default = true
        "--LSI"
            help = "Whether to use LSI or not (default = true)"
            arg_type = Bool
            required = false
            default = true
        "--max_iterations"
            help = "Number of allowed iterations before the inversions are halted (default = 10)."
            arg_type = Int
            required = false
            default = 10
        "--gamma"
            help = "Levenberg-Marquardt gamma value (default = 10)."
            arg_type = Float64
            required = false
            default = 10.0
        "--dsigma_scale"
            help = "Scale factor to (dσ)^2 to assess convergence"
            arg_type = Float64
            required = false
            default = 2.0
        "--o2_scale"
            help = "Oxygen spectroscopy scaling factor"
            arg_type = Float64
            required = false
            default = 1.0
        "--co2_scale_weak"
            help = "CO2 spectroscopy scale factor for the weak CO2 band"
            arg_type = Float64
            required = false
            default = 1.0
        "--co2_scale_strong"
            help = "CO2 spectroscopy scale factor for the strong CO2 band"
            arg_type = Float64
            required = false
            default = 1.0

    end

    return parse_args(s)

end