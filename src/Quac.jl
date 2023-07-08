module Quac

include("Gate.jl")
include("Array.jl")
include("Circuit.jl")
include("Algorithms.jl")

if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
    function __init__()
        @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("../ext/QuacMakieExt.jl")
    end
end

end
