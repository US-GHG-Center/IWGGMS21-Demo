
my_args = [
    "--solar_model",
    "./example_data/l2_solar_model.h5",
    "--L1b",
    "./example_data/2021030111564431_inputs.h5",
    "--L2Met",
    "./example_data/2021030111564431_inputs.h5",
    "--L2CPr",
    "./example_data/2021030111564431_inputs.h5",
    "--sounding_id",
    "2021030111564431",
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
    "10"
]

empty!(ARGS)

for a in my_args
    push!(ARGS, a)
end

include("run.jl")
