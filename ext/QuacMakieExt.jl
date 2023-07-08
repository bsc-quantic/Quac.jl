module QuacMakieExt

if isdefined(Base, :get_extension)
    using Quac
else
    using ..Quac
end

using Makie

end
