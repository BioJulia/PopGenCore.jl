export ishom, ishet, isbiallelic, isbinary


"""
    isbiallelic(data::GenoArray)
Returns `true` if the `GenoArray` is biallelic, `false` if not.
"""
@inline function isbiallelic(data::T) where T<:GenoArray
    length(unique(Base.Iterators.flatten(skipmissing(data)))) == 2
end

# method for using just a genodata input
function isbiallelic(data::DataFrame)
    tmp = sort(data, [:locus, :name])
    prt = Base.Iterators.partition(tmp.genotype, length(data.name.pool))
    bi = findfirst(!isbiallelic, collect(prt))
    isnothing(bi) ? true : false
end

"""
    isbiallelic(data::PopData)
Returns `true` all the loci in the `PopData` are biallelic, `false` if not.
"""
isbiallelic(data::PopData) = data.metadata.biallelic


#TODO how to treat haploids?
"""
```
ishom(locus::T) where T <: GenoArray
ishom(locus::Genotype)
ishom(locus::Missing)
```
A series of methods to test if a locus or loci are homozygous and return `true` if
it is, `false` if it isn't (or missing). For calculations, we recommend using `_ishom()`,
which returns `missing` if the genotype is `missing`. The vector version
simply maps the function over the elements.
"""
@inline function ishom(locus::Genotype)
    # if the first equals all others, return true
    return all(@inbounds first(locus) .== locus)
end

# public facing method
ishom(locus::Missing) = false
# API computational method
_ishom(locus::Missing) = missing


@inline function ishom(locus::T) where T<:GenoArray
    return @inbounds map(ishom, locus)
end

@inline function ishom(locus::T) where T<:Base.SkipMissing
    return @inbounds map(ishom, locus)
end


"""
    ishom(locus::Genotype, allele::Signed)
    ishom(loci::GenoArray, allele::Signed)
Returns `true` if the `locus`/`loci` is/are homozygous for the specified `allele`.
"""
function ishom(geno::T, allele::U) where T<:Genotype where U<:Integer
    ∈(allele, geno) & ishom(geno) ? true : false
end

ishom(geno::T, allele::U) where T<:GenoArray where U<:Integer = map(i -> ishom(i, allele), geno)

# public facing method
ishom(geno::Missing, allele::U) where U<:Integer = false
# API computational method
_ishom(geno::Missing, allele::U) where U<:Integer = missing

"""
```
ishet(locus::T) where T <: GenoArray
ishet(locus::Genotype)
ishet(locus::Missing)
```
A series of methods to test if a locus or loci are heterozygous and return `true` if
it is, `false` if it isn't. The vector version simply broadcasts the function over the
elements. Under the hood, this function is simply `!ishom`.
"""
@inline function ishet(locus::Genotype)
    return !ishom(locus)
end

# public facing method
ishet(locus::Missing) = false
# API computational method
_ishet(locus::Missing) = missing


@inline function ishet(locus::T) where T<:GenoArray
    return @inbounds map(ishet, locus)
end


@inline function ishet(locus::T) where T<:Base.SkipMissing
    return @inbounds map(ishet, locus)
end


"""
    ishet(locus::Genotype, allele::Signed)
    ishet(loci::GenoArray, allele::Signed)
Returns `true` if the `locus`/`loci` is/are heterozygous for the specified `allele`. 
"""
function ishet(geno::T, allele::U) where T<:Genotype where U<:Integer
    ∈(allele, geno) & !ishom(geno) ? true : false
end

ishet(geno::T, allele::U) where T<:GenoArray where U<:Integer = map(i -> ishet(i, allele), geno)
# public facing method
ishet(geno::Missing, allele::U) where U<:Integer = false
# API computational method
_ishet(geno::Missing, allele::U) where U<:Integer = missing

