#=

    Assorted functions to help with sending and receiving data via MPI

=#


# Which are the variable names we want to send/recieve via MPI?
MPI_variables_to_send = [
    :scene_inputs
]


# Generate a tag dictionary for each symbol, accessed by the MPI rank number:
# e.g. MPI_tags[:scene_inputs] = Dict(1 => 100)
# means: variable "scene_inputs" will be communicated to rank 1 with tag 100

MPI_tags = Dict{Symbol, Dict{Int, Int}}()

for (cnt,sym) in enumerate(MPI_variables_to_send)

    MPI_tags[sym] = Dict{Int,Int}()

    start_tag = (100 + (cnt-1) * (MPI_size * length(MPI_variables_to_send)))
    for rank in 1:MPI_size-1
        MPI_tags[sym][rank] = start_tag + rank - 1
    end
end