"""
    isallmissing(itr::T) where T<:AbstractArray
    isallmissing(itr::T) where T<:Tuple
A lightning-fast and non-allocating conditional which tests if iterable `itr`
contains only `missing` values.
"""
function isallmissing(itr::T) where T <: AbstractArray
    val = true
    i = 1
    while (i <= length(itr)) & val
        val = itr[i] === missing
        i += 1
    end
    return val
end

function isallmissing(itr::T) where T <: Tuple
    val = true
    i = 1
    while (i <= length(itr)) & val
        val = itr[i] === missing
        i += 1
    end
    return val
end


"""
    nonmissing(vec::T) where T<:AbstractArray
Convenience function to count the number of non-`missing` values
in a vector.
"""
@inline function nonmissing(vec::T) where T<:AbstractArray
    mapreduce(!ismissing, +, vec)
end


"""
    nonmissing(data::PopData, locus::String)
Convenience function to count the number of non-`missing` samples
at a locus.
"""
@inline function nonmissing(data::PopData, locus::String)
    data.genodata[data.genodata.locus .== locus, :genotype] |> nonmissing
end


"""
    nonmissings(vec1::AbstractVector, vec2::AbstractVector)
Return a vector of indices where neither input vectors have a `missing` value, i.e. an
intersection of the indices of their non-missing elements.
"""
@inline function nonmissings(vec1::T, vec2::T) where T <: AbstractVector
    mapreduce(i -> findall(!ismissing, i), intersect, (vec1, vec2))
end