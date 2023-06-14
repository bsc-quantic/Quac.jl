module Quac

using Requires: @require

include("Gate.jl")
include("Array.jl")
include("Circuit.jl")
include("Algorithms.jl")
include("parser.jl")

function __init__()
    @require Luxor = "ae8d54c2-7ccd-5906-9d76-62fc9837b5bc" include("Visualization.jl")
    @require VSCodeServer = "9f5989ce-84fe-42d4-91ec-6a7a8d53ed0f" using Luxor
end

end
