"""
    isbiallelic(data::GenoArray)
Returns `true` if the `GenoArray` is biallelic, `false` if not.
"""
@inline function isbiallelic(data::T) where T<:GenoArray
    length(unique(Base.Iterators.flatten(skipmissing(data)))) == 2
end
#precompile(isbiallelic, (Vector{SNP}))
#precompile(isbiallelic, (Vector{MSat}))

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
ishom(locus::T) where T<:Base.SkipMissing
ishom(locus::Genotype)
ishom(locus::Missing)
```
A series of methods to test if a locus or loci are homozygous and return `true` if
it is, `false` if it isn't (or missing). For calculations, we recommend using `_ishom()`,
which returns `missing` if the genotype is `missing`. The vector methods
simply map the function over the elements. Haploid genotypes return `false`.
"""
ishom(locus::Tuple{Int8})::Bool = false 
ishom(locus::Tuple{Int16})::Bool = false

function ishom(locus::NTuple{N,Int8})::Bool where N
    out = true
    @inbounds a1 = locus[1]
    i = 2
    while i <= N
        @inbounds out = a1 == locus[i]
        !out && break
        i += 1
    end
    return out
end

function ishom(locus::NTuple{N,Int16})::Bool where N
    out = true
    @inbounds a1 = locus[1]
    i = 2
    while i <= N
        @inbounds out = a1 == locus[i]
        !out && break
        i += 1
    end
    return out
end

# public facing method
ishom(locus::Missing)::Bool = false
ishom(locus::T)::Vector{Bool} where T<:GenoArray = @inbounds map(ishom, locus)
ishom(locus::T)::Vector{Bool} where T<:Base.SkipMissing = @inbounds map(ishom, locus)
precompile(ishom, (NTuple{2,Int8},))
precompile(ishom, (NTuple{2,Int16},))
precompile(ishom, (Vector{Union{Missing,SNP}},))
precompile(ishom, (Vector{Union{Missing,MSat}},))


# API computational methods
_ishom(locus::NTuple{N,Int8})::Bool where N = ishom(locus)

_ishom(locus::NTuple{N,Int16})::Bool where N = ishom(locus)

_ishom(locus::T)::Vector{Bool} where T<:GenoArray = @inbounds map(_ishom, locus)
_ishom(locus::T)::Vector{Bool} where T<:Base.SkipMissing = @inbounds map(_ishom, locus)
_ishom(locus::Missing)::Missing = missing
precompile(_ishom, (NTuple{2,Int8},))
precompile(_ishom, (NTuple{2,Int16},))
precompile(_ishom, (Vector{Union{Missing,SNP}},))
precompile(_ishom, (Vector{Union{Missing,MSat}},))

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
ishom(geno::Missing, allele::U)::Bool where U<:Integer = false

# API computational method
function _ishom(geno::T, allele::U)::Bool where T<:Genotype where U<:Integer
    !_ishom(geno) && return false 
    ∈(allele, geno) && return true
end
function _ishom(geno::T, allele::U)::Vector{Bool} where T<:GenoArray where U<:Integer
    @inbounds [ishom(i, allele) for i in geno]
end
_ishom(geno::Missing, allele::U)::Bool where U<:Integer = missing
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
ishet(locus::T)::Vector{Bool} where T<:GenoArray = @inbounds map(ishet, locus)
ishet(locus::T)::Vector{Bool} where T<:Base.SkipMissing = @inbounds map(ishet, locus)

_ishet(locus::Genotype) = !_ishom(locus)
_ishet(locus::Missing) = missing
_ishet(locus::T)::Vector{Bool} where T<:GenoArray = @inbounds map(_ishet, locus)
_ishet(locus::T)::Vector{Bool} where T<:Base.SkipMissing = @inbounds map(_ishet, locus)


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

ishet(geno::Missing, allele::U)::Bool where U<:Integer = false

function _ishet(geno::T, allele::U)::Bool where T<:Genotype where U<:Integer
    _ishom(geno) && return false
    ∈(allele, geno) && return true
end

function _ishet(geno::T, allele::U)::Vector{Bool} where T<:GenoArray where U<:Integer
    [ishet(i, allele) for i in geno]
end

_ishet(geno::Missing, allele::U)::Bool where U<:Integer = missing


