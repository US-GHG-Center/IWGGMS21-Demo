using Revise # May want to use this script for debugging..

snid = "2021061809241801"
#snid = "2021061809234403"


my_args = [
    "--solar_model",
    "./data/hybrid_reference_spectrum_p005nm_resolution_c2022-11-30_with_unc.nc",
    #"data/TSIS-SIM-scaled_l2_solar_model.h5",
    "--L1b",
    "./example_data/2021030111564431_inputs.h5",
    #"../ACOS-validation/0_Belchatow_Nassar_1/oco3_L1bScSC_12023a_210618_B11064r_240526223245.h5",
    #"data/oco2_L1bScTG_35448a_210301_B11006r_220330101215.h5",
    #"data/oco2_L1bScND_42294a_220614_B11008r_220810203728.h5",
    "--L2Met",
    "./example_data/2021030111564431_inputs.h5",
    #"../ACOS-validation/0_Belchatow_Nassar_1/oco3_L2MetSC_12023a_210618_B11064r_240526223147.h5",
    #"data/oco2_L2MetTG_35448a_210301_B11006r_220327161030.h5",
    #"data/oco2_L2MetND_42294a_220614_B11008r_220810210117.h5",
    "--L2CPr",
    "./example_data/2021030111564431_inputs.h5",
    #"../ACOS-validation/0_Belchatow_Nassar_1/oco3_L2CPrSC_12023a_210618_B11064r_240526224849.h5",
    #"data/oco2_L2CPrTG_35448a_210301_B11006r_220327173727.h5",
    #"data/oco2_L2CPrND_42294a_220614_B11008r_220811060138.h5",
    "--sounding_id",
    #snid,
    #"2021061809234403",
    "2021030111564431",
    #"2021030112001905",
    #"2022061414300438",
    "--spec",
    "1,2,3",
    "--aerosols",
    "true",
    "--o2_scale",
    "1.0048",
    "--co2_scale_weak",
    "0.994",
    "--co2_scale_strong",
    "0.998",
    "--gamma",
    "1000.0",
    "--dsigma_scale",
    "1.0",
    "--max_iterations",
    "10",
    "--output",
    "$(snid).h5",
    #"2021030111564431_two_os.h5",
    #"2021030112001905.h5",
    #"2022061414300438.h5"
]

empty!(ARGS)

for a in my_args
    push!(ARGS, a)
end

include("run.jl")
