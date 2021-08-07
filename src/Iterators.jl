export skipinf, skipnan, skipinfnan

"""
    nonmissings(vec1::AbstractVector, vec2::AbstractVector)
Return a vector of indices where neither input vectors have a `missing` value, i.e. an
intersection of the indices of their non-missing elements.
"""
@inline function nonmissings(vec1::T, vec2::T) where T <: AbstractVector
    mapreduce(i -> findall(!ismissing, i), intersect, (vec1, vec2))
end


"""
skipnan(itr)
Return an iterator over the elements in `itr` skipping `Inf` and `-Inf` values. The returned
object can be indexed using indices of itr if the latter is indexable. Indices
corresponding to `Inf` values are not valid: they are skipped by keys and eachindex,   
and a MissingException is thrown when trying to use them. This is effectively `skipmissing`
for `Inf` and `-Inf` values.

Use collect to obtain an `Array` containing the non-`Inf` values in `itr`. Note that even  
if `itr` is a multidimensional array, the result will always be a `Vector` since it is not   
possible to remove `Inf`s while preserving dimensions of the input.
"""
skipinf(itr) = Iterators.filter(isfinite, itr)


"""
skipnan(itr)
Return an iterator over the elements in `itr` skipping `NaN` values. The returned
object can be indexed using indices of itr if the latter is indexable. Indices
corresponding to `NaN` values are not valid: they are skipped by keys and eachindex,   
and a MissingException is thrown when trying to use them. This is effectively `skipmissing`
for `NaN` values.

Use collect to obtain an `Array` containing the non-`NaN` values in `itr`. Note that even  
if `itr` is a multidimensional array, the result will always be a `Vector` since it is not   
possible to remove `NaN`s while preserving dimensions of the input.
"""
skipnan(itr) = Iterators.filter(!isnan, itr)

"""
skipinfnan(itr)
Return an iterator over the elements in `itr` skipping `NaN`, `Inf` and `-Inf` values.
See the docstrings of `skipinf` and `skipnan` more details.
"""
skipinfnan(itr) = Iterators.filter(x -> (isfinite(x) & !isnan(x)), itr)