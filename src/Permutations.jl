export permute_loci!, permute_samples!, permute_genotypes!, permute_alleles!
export strict_shuffle, strict_shuffle!

"""
    permute_loci!(data::PopData)
Edits `PopData` in place with loci permuted across populations within
the `.loci` dataframe.
"""
@inline function permute_loci!(data::PopData)
    @inbounds @sync for locus in groupby(data.genodata, :locus)
        Base.Threads.@spawn begin
            shuffle!(locus.population)
        end
    end
    data
end

"""
    permute_samples!(data::PopData; meta::Bool = false)
Edits `PopData` in place with samples permuted across populations within
the `.loci` dataframe. Since performance is important for many permutations,
the default is to only edit the `.loci` table in place; use `meta = true`
if you also require the `.meta` dataframe edited in place.
"""
@inline function permute_samples!(data::PopData; meta::Bool = false)
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


@inline function permute_samples!(data::AbstractDataFrame, popnames::Vector{String})
    pops = shuffle(popnames)

    @inbounds #=@sync=# for name in groupby(data, :name)
        #Base.Threads.@spawn begin 
            @inbounds name.population .= pop!(pops)
        #end
    end
    data
end


"""
    permute_genotypes!(data::PopData; by::String = "locus")
Edits `PopData` in place with genotypes permuted across individuals within
the `.loci` dataframe. Use `by = "population"` to permute genotypes
within populations.
"""
@inline function permute_genotypes!(data::PopData; by::String = "locus")
    # establish mode of operation
    groupings = by == "locus" ? :locus : [:locus, :population]
    @inbounds @sync for grp in groupby(data.genodata, groupings)
        Base.Threads.@spawn begin
            grp.genotype .= strict_shuffle!(grp.genotype)
        end
    end
    data
end


"""
    permute_alleles!(data::PopData; ploidy::Union{Nothing, Int} = nothing, by::String = "locus")
Edits `PopData` in place with alleles permuted and reconstructed into genotypes
for each locus within the `.loci` dataframe. Use `by = "population"`
to permute alleles within populations. If `ploidy` is not provided (default `ploidy = nothing`),
then ploidy will be identified from the PopData. If performance is important,
it would be best to identify ploidy in advance and set it to a specific integer.
"""
@inline function permute_alleles!(data::PopData; ploidy::Union{Nothing, Int} = nothing, by::String = "locus")
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
    strict_shuffle(x::T) where T <: AbstractArray
Shuffle only the non-missing values of a Vector and return a copy of the vector,
keeping the `missing` values at their original locations.
Use `strict_shuffle!` to edit in-place instead of returning a copy.
"""
@inline function strict_shuffle(x::T) where T <: AbstractArray
    y = copy(x)
    @inbounds shuffle!(@view y[.!ismissing.(y)])
    return y
end

"""
    strict_shuffle!(x::T)! where T <: AbstractArray
Shuffle only the non-missing values of a Vector, keeping the
`missing` values at their original locations. Use `strict_shuffle`
to return a copy instead of editing in-place.
"""
@inline function strict_shuffle!(x::T) where T <: AbstractArray
    @inbounds shuffle!(@view x[.!ismissing.(x)])
    return x
end