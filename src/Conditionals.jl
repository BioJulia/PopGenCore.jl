"""
    isbiallelic(data::GenoArray)
Returns `true` if the `GenoArray` is biallelic, `false` if not.
"""
@inline function isbiallelic(data::T) where T<:GenoArray
    allelecount(data) == 2
end
#precompile(isbiallelic, (Vector{SNP}))
#precompile(isbiallelic, (Vector{MSat}))

# method for using just a genodata input
function isbiallelic(data::DataFrame)
    tmp = sort(data, :locus)
    @inbounds for i in Base.Iterators.partition(tmp.genotype, length(data.name.pool))
        isbiallelic(i) || return false
    end
    return true
end

"""
    isbiallelic(data::PopData)
Returns `true` all the loci in the `PopData` are biallelic, `false` if not.
"""
isbiallelic(data::PopData) = data.metadata.biallelic


# public facing method
# Haploids treated as heterozygotes
"""
```
ishom(locus::T) where T <: GenoArray
ishom(locus::T) where T<:Base.SkipMissing
ishom(locus::Genotype)
ishom(locus::Missing)
```
A series of methods to test if a locus or loci are homozygous and return `true` if
it is, `false` if it isn't (or missing). For calculations, we recommend using `_ishom()`,
which returns `missing` if the genotype is `missing`. The vector methods
simply map the function over the elements. Haploid genotypes return `false`.
"""
function ishom(geno::NTuple{N,T})::Bool where N where T<:Union{Int8, Int16}
    first,rest = Base.Iterators.peel(geno)
    @inbounds for v in rest
        v != first && return false
    end
    return true
end

ishom(geno::NTuple{1,T}) where T<:Union{Int8, Int16} = false

ishom(locus::Missing)::Bool = false
function ishom(locus::T)::Vector{Bool} where T<:GenoArray
    @inbounds map(ishom, locus)
end
function ishom(locus::T)::Vector{Bool} where T<:Base.SkipMissing
    @inbounds map(ishom, locus)
end
precompile(ishom, (NTuple{2,Int8},))
precompile(ishom, (NTuple{2,Int16},))
precompile(ishom, (Vector{Union{Missing,SNP}},))
precompile(ishom, (Vector{Union{Missing,MSat}},))


# API computational methods
function _ishom(geno::NTuple{N,T})::Bool where N where T<:Union{Int8, Int16}
    ishom(geno)
end
_ishom(locus::Missing) = missing
_ishom(locus::T) where T<:GenoArray = @inbounds map(_ishom, locus)
_ishom(locus::T) where T<:Base.SkipMissing = @inbounds map(_ishom, locus)

#= scales for size, which isn't super necessary yet
=#
"""
    ishom(locus::Genotype, allele::Integer)
    ishom(loci::GenoArray, allele::Integer)
Returns `true` if the `locus`/`loci` is/are homozygous for the specified `allele`.
"""
function ishom(geno::T, allele::U)::Bool where T<:Genotype where U<:Integer
    !ishom(geno) && return false
    ∈(allele, geno) && return true
end
precompile(ishom, (NTuple{2,Int8},Int64))
precompile(ishom, (NTuple{2,Int8},Int8))
precompile(ishom, (NTuple{2,Int16},Int64))
precompile(ishom, (NTuple{2,Int16},Int16))

function ishom(geno::T, allele::U)::Vector{Bool} where T<:GenoArray where U<:Integer 
    @inbounds [ishom(i, allele) for i in geno]
end

# public facing method
ishom(geno::Missing, allele::U) where U<:Integer = false

# API computational method
function _ishom(geno::T, allele::U)::Bool where T<:Genotype where U<:Integer
    !_ishom(geno) && return false 
    ∈(allele, geno) && return true
end
function _ishom(geno::T, allele::U) where T<:GenoArray where U<:Integer
    @inbounds [ishom(i, allele) for i in geno]
end
_ishom(geno::Missing, allele::U) where U<:Integer = missing
precompile(_ishom, (NTuple{2,Int8},Int8))
precompile(_ishom, (NTuple{2,Int8},Int64))
precompile(_ishom, (NTuple{2,Int16},Int64))
precompile(_ishom, (NTuple{2,Int16},Int16))

# public facing method
"""
```
ishet(locus::T) where T <: GenoArray
ishet(locus::Genotype)
ishet(locus::Missing)
```
A series of methods to test if a locus or loci are heterozygous and return `true` if
it is, `false` if it isn't (or missing). For calculations, we recommend using `_ishet()`,
which returns `missing` if the genotype is `missing`. The vector version
simply maps the function over the elements.
"""
ishet(locus::Genotype) = !ishom(locus)
ishet(locus::Missing) = false
function ishet(locus::T)::Vector{Bool} where T<:GenoArray
    @inbounds map(ishet, locus)
end
function ishet(locus::T)::Vector{Bool} where T<:Base.SkipMissing
    @inbounds map(ishet, locus)
end
_ishet(locus::Genotype) = !_ishom(locus)
_ishet(locus::Missing) = missing
function _ishet(locus::T) where T<:GenoArray
    @inbounds map(_ishet, locus)
end
function _ishet(locus::T) where T<:Base.SkipMissing
    @inbounds map(_ishet, locus)
end

# public facing methods
"""
    ishet(locus::Genotype, allele::Signed)
    ishet(loci::GenoArray, allele::Signed)
Returns `true` if the `locus`/`loci` is/are heterozygous for the specified `allele`. 
"""
function ishet(geno::T, allele::U)::Bool where T<:Genotype where U<:Integer
    ishom(geno) && return false
    ∈(allele, geno) && return true
end

function ishet(geno::T, allele::U)::Vector{Bool} where T<:GenoArray where U<:Integer
    [ishet(i, allele) for i in geno]
end

ishet(geno::Missing, allele::U) where U<:Integer = false

function _ishet(geno::T, allele::U)::Bool where T<:Genotype where U<:Integer
    _ishom(geno) && return false
    ∈(allele, geno) && return true
end

function _ishet(geno::T, allele::U) where T<:GenoArray where U<:Integer
    [ishet(i, allele) for i in geno]
end

_ishet(geno::Missing, allele::U) where U<:Integer = missing


