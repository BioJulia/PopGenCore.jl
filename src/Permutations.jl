"""
    permuteloci!(data::PopData)
Edits `PopData` in place with loci permuted across populations within
the `.genodata` dataframe.
"""
@inline function permuteloci!(data::PopData)
    @inbounds @sync for locus in groupby(data.genodata, :locus)
        Base.Threads.@spawn begin
            shuffle!(locus.population)
        end
    end
    data
end

"""
    permutesamples!(data::PopData; meta::Bool = false)
Edits `PopData` in place with samples permuted across populations within
the `.genodata` dataframe. Since performance is important for many permutations,
the default is to only edit the `.genodata` table in place; use `meta = true`
if you also require the `.meta` dataframe edited in place.
"""
@inline function permutesamples!(data::PopData; meta::Bool = false)
    pops = shuffle(data.sampleinfo.population)

    if meta == true
        meta_pops = copy(pops)
        @inbounds for name in groupby(data.metadata, :name)
            @inbounds name.population .= pop!(meta_pops)
        end
    end
    @inbounds @sync for name in groupby(data.genodata, :name)
        Base.Threads.@spawn begin
            @inbounds name.population .= pop!(pops) 
        end
    end
    data
end


@inline function permutesamples!(data::AbstractDataFrame, popnames::Vector{String})
    pops = shuffle(popnames)

    @inbounds #=@sync=# for name in groupby(data, :name)
        #Base.Threads.@spawn begin 
            @inbounds name.population .= pop!(pops)
        #end
    end
    data
end


"""
    permutegenotypes!(data::PopData; by::String = "locus")
Edits `PopData` in place with genotypes permuted across individuals within
the `.genodata` dataframe. Use `by = "population"` to permute genotypes
within populations.
"""
@inline function permutegenotypes!(data::PopData; by::String = "locus")
    # establish mode of operation
    groupings = by == "locus" ? :locus : [:locus, :population]
    @inbounds @sync for grp in groupby(data.genodata, groupings)
        Base.Threads.@spawn begin
            grp.genotype .= strictshuffle!(grp.genotype)
        end
    end
    data
end


"""
    permutealleles!(data::PopData; ploidy::Union{Nothing, Int} = nothing, by::String = "locus")
Edits `PopData` in place with alleles permuted and reconstructed into genotypes
for each locus within the `.genodata` dataframe. Use `by = "population"`
to permute alleles within populations. If `ploidy` is not provided (default `ploidy = nothing`),
then ploidy will be identified from the PopData. If performance is important,
it would be best to identify ploidy in advance and set it to a specific integer.
"""
@inline function permutealleles!(data::PopData; ploidy::Union{Nothing, Int} = nothing, by::String = "locus")
    if ploidy === nothing
        length(data.metadata.ploidy) > 1 && error("This permutation method is not appropriate for mixed-ploidy data")
        ploidy = first(tmp)
    end

    # establish mode of operation
    groupings = by == "locus" ? :locus : [:locus, :population]

    @inbounds @sync for grp in groupby(data.genodata, groupings)
        Base.Threads.@spawn begin
            alle = shuffle(alleles(grp.genotype))
            new_genos = Tuple.(Base.Iterators.partition(alle, ploidy))
            (@view grp.genotype[.!ismissing.(grp.genotype)]) .= new_genos
        end
    end
    data
end


"""
    strictshuffle(x::T) where T <: AbstractArray
Shuffle only the non-missing values of a Vector and return a copy of the vector,
keeping the `missing` values at their original locations.
Use `strictshuffle!` to edit in-place instead of returning a copy.
"""
@inline function strictshuffle(x::T) where T <: AbstractArray
    y = copy(x)
    @inbounds shuffle!(@view y[.!ismissing.(y)])
    return y
end

"""
    strictshuffle!(x::T)! where T <: AbstractArray
Shuffle only the non-missing values of a Vector, keeping the
`missing` values at their original locations. Use `strictshuffle`
to return a copy instead of editing in-place.
"""
@inline function strictshuffle!(x::T) where T <: AbstractArray
    @inbounds shuffle!(@view x[.!ismissing.(x)])
    return x
end