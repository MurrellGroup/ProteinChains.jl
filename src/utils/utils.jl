include("renumber.jl")
export renumber!

include("truncate.jl")
export truncate

include("repeats.jl")
export detect_repeats

include("secondary.jl")

method(mmcifdict::BioStructures.MMCIFDict) = get(mmcifdict, "_exptl.method", ["missing"]) |> first
