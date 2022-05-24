"""
    partitionarray(array::AbstractArray, steps::AbstractVector{<:Integer})
Like Base.Iterators.Partition, except you can apply arbitrary sizes to
partition the array by. The `steps` must add up to the total row length
of the array.

**Example**
```
julia> partitionarray(rand(20,5), [10,3,4,3]) .|> size
[(10, 5), (3, 5), (4, 5), (3, 5)]
```
"""
function partitionarray(array::AbstractArray, steps::AbstractVector{<:Integer})
    # solution brilliantly provided by @stevengj and @mcabbott on Slack and Discourse 
    # (https://discourse.julialang.org/t/is-there-a-simple-intuitive-way-to-partition-a-matrix-by-arbitrary-strides-like-i/55863)
    v = axes(array,1)
    v == 1:sum(steps) || error("Steps provided do not sum to length of the first dimension")
    i = firstindex(v)
    tmp = (view(v, i:(i+=s)-1) for s in steps)
    [view(array,r,:) for r in tmp]
end


"""
    pairwisepairs(smp_names::Vector{T}) where T
Given a vector, returns a lazy iterator of tuples of unique all x 
all combinations of element pairs, excluding self-comparisons.

**Example**
```
julia> colors = ["red_1", "red_2", "blue_1", "blue_2"] ;

julia> pairwisepairs(colors) |> collect
6-element Array{Tuple{String,String},1}:
 ("red_1", "red_2")
 ("red_1", "blue_1")
 ("red_1", "blue_2")
 ("red_2", "blue_1")
 ("red_2", "blue_2")
 ("blue_1", "blue_2")
```
"""
@inline function pairwisepairs(smp_names::AbstractVector{T}) where T
    len = length(smp_names)
    (tuple(smp_names[i], smp_names[j]) for i in 1:len-1 for j in i+1:len)
end

"""
    simpairs(data::Vector{String})
Takes a Vector of sample names and returns a Tuple of sample pairs, grouped by simulation
number. This is an internal function used for isolating sibship pairs from simulated shipship
pairs (via `PopGenSims.jl`) to perform `relatedness` estimates only on those pairs.

**Example**
```julia
julia> a = ["sim1_off1", "sim1_off2", "sim2_off1", "sim2_off2"] ;

julia> _simpairs(a)
("sim1_off1", "sim1_off2")
("sim2_off1", "sim2_off2")
```
"""
function simpairs(data::Vector{String})
    n = length(data)
    isodd(n) && throw(ArgumentError("Expected an even number of samples, but got $n"))
    Tuple.(Base.Iterators.partition(sort(data), 2))
end


"""
    skipinf(itr)
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
