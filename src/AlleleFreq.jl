export allele_freq, allele_freq_vec, avg_allele_freq

"""
    allele_freq(data::PopData)
Return a NamedTuple of `Dicts` of allele frequencies of all
loci in a `PopData` object. These are global allele frequencies.
"""
function allele_freq(data::PopData)
    NamedTuple{Tuple(Symbol.(loci(data)))}(
            DataFrames.combine(
            groupby(data.genodata, :locus),
            :genotype => allele_freq => :frq
        )[:, :frq] |> Tuple
    )
end

"""
    allele_freq(locus::GenoArray)
Return a `Dict` of allele frequencies of a GenoArray (typically a single locus) in 
a `PopData` object.
"""
@inline function allele_freq(locus::GenoArray)
    all(ismissing.(locus)) == true && return Dict{eltype(nonmissingtype(eltype(locus))), Float64}()
    proportionmap(alleles(locus))
end

"""
    allele_freq(geno::Genotype)
Return a `Dict` of allele frequencies of the alleles within a single Genotype in a `PopData`
object.
"""
@inline function allele_freq(geno::Genotype)
    d = Dict{eltype(geno),Float32}()
    len = length(geno)
    @inbounds @simd for allele in geno
        d[allele] = @inbounds get!(d, allele, 0.0) + 1.0/len
    end
    return d
end

"""
    allele_freq_vec(locus::GenoArray)
Return a Vector of allele frequencies of a single locus in a `PopData`
object. Similar to `allele_freq()`, except it returns only the frequencies,
without the allele names, meaning they can be in any order. This is useful
for getting the expected genotype frequencies.
"""
@inline function allele_freq_vec(locus::GenoArray)
    flat_alleles = alleles(locus)
    len = length(flat_alleles)
    d = [count(i -> i == j, flat_alleles) for j in unique(flat_alleles)]
    return d ./ len
end

@inline function allele_freq_vec(::Missing)
    return missing
end


"""
    avg_allele_freq(allele_dicts::AbstractVector{Dict{T, Float64}}, power::Int = 1)
Takes a Vector of Dicts generated by `allele_freq` and returns a Dict of the average
allele frequencies raised to the `power` (exponent) specified (default: `1`). 
This is typically done to calculate average allele frequencies across populations.

**Example**
```
cats = @nancycats;
alleles_df = DataFrames.combine(
    groupby(cats.loci, [:locus, :population]),
    :genotype => allele_freq => :alleles
);
DataFrames.combine(
    groupby(alleles_df, :locus),
    :alleles => (i -> sum(avg_allele_freq(i, 2))) => :avg_freq
)
```
"""
function avg_allele_freq(allele_dicts::AbstractVector{Dict{T, Float64}}, power::Int = 1) where T<:Integer   
   sum_dict = Dict{Int16, Tuple{Float32, Int}}()
   # remove any dicts with no entries (i.e. from a group without that locus)
   allele_dicts = allele_dicts[findall(i -> length(i) > 0, allele_dicts)]
   # create a list of all the alleles
   all_alleles = keys.(allele_dicts) |> Base.Iterators.flatten |> collect |> unique
   # populate the sum dict with allele frequency and n for each allele
   @inbounds for allele in all_alleles
       for allele_dict in allele_dicts
           @inbounds sum_dict[allele] = get!(sum_dict, allele, (0., 0)) .+ (get!(allele_dict, allele, 0.), 1)
       end
   end
   avg_dict = Dict{Int16, Float32}()
   # collapse the sum dict into a dict of averages
   @inbounds for (key, value) in sum_dict
       freq_sum, n = value
       avg = (freq_sum / n) ^ power
       # drop zeroes
       if !iszero(avg)
           @inbounds avg_dict[key] = avg
       end
   end
   return avg_dict
end

# method for nei_fst (pairwise)
function avg_allele_freq(allele_dicts::T, power::Int = 1) where T<:Tuple
    sum_dict = Dict{Int16, Tuple{Float32, Int}}()
    # remove any dicts with no entries (i.e. from a group without that locus)
    allele_dicts = allele_dicts[findall(i -> length(i) > 0, allele_dicts)]
    # create a list of all the alleles
    all_alleles = keys.(allele_dicts) |> Base.Iterators.flatten |> collect |> unique
    # populate the sum dict with allele frequency and n for each allele
    @inbounds for allele in all_alleles
        for allele_dict in allele_dicts
            @inbounds sum_dict[allele] = get!(sum_dict, allele, (0., 0)) .+ (get!(allele_dict, allele, 0.), 1)
        end
    end
    avg_dict = Dict{Int16, Float32}()
    # collapse the sum dict into a dict of averages
    @inbounds for (key, value) in sum_dict
        freq_sum, n = value
        avg = (freq_sum / n) ^ power
        # drop zeroes
        if !iszero(avg)
            @inbounds avg_dict[key] = avg
        end
    end
    return avg_dict
 end


"""
    allele_freq(data::PopData, locus::String; population::Bool = false)
Return a `Dict` of allele frequencies of a single locus in a `PopData`
object. Use `population = true` to return a table of allele frequencies
of that locus per population.
### Example
```
cats = @nancycats;
allele_freq(cats, "fca8")
allele_freq(cats, "fca8", population = true)
```
"""
function allele_freq(data::PopData, locus::String; population::Bool=false)
    if !population
        data.genodata[data.genodata.locus .== locus, :genotype] |> allele_freq
    else
        tmp = groupby(data.genodata[data.genodata.locus .== locus, :], :population)
        DataFrames.combine(tmp, :genotype => allele_freq => :frequency)
    end
end


#TODO swtich order of args do it's data, allele?
# Does doing that break anything?
"""
    allele_freq(allele::Int, genos::GenoArray)
Return the frequency of a specific `allele` from a vector of `genotypes`.

### Example
```
using DataFramesMeta
ncats = @nancycats;
ncats_sub @where(ncats.loci, :locus .== "fca8", :genotype .!== missing)
pop_grp = groupby(ncats_sub, :population)
DataFrames.combine(pop_grp, :genotype => (geno -> allele_freq(137, geno)) => :freq_137)
```
"""
function allele_freq(allele::Int, genos::T) where T<:GenoArray
    tmp = allele_freq(genos)
    haskey(tmp, allele) ? getindex(tmp, allele) : 0.
end