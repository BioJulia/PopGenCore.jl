export reciprocal, reciprocal_sum

## Utilities relevant for math ##

"""
reciprocal(num::T) where T <: Signed
Returns the reciprocal (1/number) of a number. Will return `0` when
the number is `0` instead of returning `Inf`.
"""
function reciprocal(num::T) where T <: Real
    !iszero(num) ? 1.0/float(num) : 0.0
end


"""
reciprocal_sum(x::AbstractVector{T}) where T<:Real
Return the sum of the reciprocal values of `x`, skipping the `Inf` values
resulting from divide-by-zero errors.
"""
function reciprocal_sum(x::AbstractVector{T}) where T<:Real
    mapreduce(reciprocal, +, x)
end