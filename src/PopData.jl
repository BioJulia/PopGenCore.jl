"""
    AbstractType PopObj
Generic AbstractType for use in PopGen.jl
"""
abstract type PopObj end


"""
```
mutable struct PopDataInfo
    samples::Int64
    sampeinfo::DataFrame
    loci::Int64
    locusinfo::DataFrame
    populations::Int64
    ploidy::Union{Int8, Vector{Int8}}
    biallelic::Bool
end
```
The data struct used internally as `PopData.metadata` fields to store basic information
about the `PopData` for easy access.
"""
mutable struct PopDataInfo
    samples::Int64
    sampleinfo::DataFrame
    loci::Int64
    locusinfo::DataFrame
    populations::Int64
    ploidy::Union{Int8, Vector{Int8}}
    biallelic::Bool
end

# constructor FORMAT just the genodata dataframe
function PopDataInfo(genodf::DataFrame)
    nameorder = unique(genodf.name)
    sampleinfo = unique(dropmissing(genodf), :name)
    # restore correct sample order if first record for a sample has missing genotype
    sampleinfo = sampleinfo[indexin(nameorder,sampleinfo.name),:]
    sampleinfo.ploidy = [ismissing(geno) ? Int8(0) : Int8(length(geno)) for geno in sampleinfo.genotype]
    select!(sampleinfo, :name => collect => :name, :population, :ploidy)
    ploidy = unique(sampleinfo.ploidy)
    ploidy = length(ploidy) == 1 ? ploidy[1] : ploidy
    PopDataInfo(
        length(genodf.name.pool),
        sampleinfo,
        length(genodf.locus.pool),
        DataFrame(:chromosome => PooledArray(fill("0",length(genodf.locus.pool)) , compress = true), :locus => genodf.locus.pool, :cm => Int8(0), :bp => Int8(0)),
        length(genodf.population.pool),
        ploidy,
        isbiallelic(genodf)
    )
end

function Base.show(io::IO, data::PopDataInfo)
    println(io, " ploidy:       ", join(data.ploidy, ","))
    println(io, " loci:         ", data.loci)
    println(io, " samples:      ", data.samples)
    println(io, " populations:  ", data.populations)
    println(io, " biallelic:    ", data.biallelic)
    println(io, " other fields: locusinfo, sampleinfo")
end

"""
```
PopData
    metadata::PopDataInfo
    genodata::DataFrame
```
The data struct used for the PopGen population genetics ecosystem. You are
**strongly** discouraged from manually creating tables to pass into a PopData,
and instead should use the provided file importers and utilities.

    - `metadata` PopDataInfo of  data information
        - `samples` - the number of samples in the data
        - `sampleinfo` - DataFrame of sample names,populations, ploidy, etc.
        - `loci` - the number of loci in the data
        - `locusinfo` - DataFrame of locus names, chromosome, physical position, etc.
        - `populations` - the number of populations in the data
        - `ploidy` - the ploidy (or ploidies) present in the data
        - `biallelic` - if all the markers are biallelic
    - `genodata` DataFrame of sample genotype records
        - `name` - the individual/sample names [`PooledArray`]
        - `population` - population names [`PooledArray`]
        - `locus` - locus names [`PooledArray`]
        - `genotype` - genotype values [`NTuple{N,Signed}`]

"""
struct PopData <: PopObj
    metadata::PopDataInfo
    genodata::DataFrame
    function PopData(meta::PopDataInfo, loci::DataFrame)
        if "ploidy" ∉ names(meta.sampleinfo)
            insertcols!(meta.sampleinfo, 3, :ploidy => Int8(0))
        end
        metacols = ["name", "population", "ploidy"]
        metacheck = intersect(metacols, names(meta.sampleinfo))
        if metacheck != metacols
            throw(error("metadata sampleinfo missing columns $(symdiff(metacheck, metacols))"))
        end
        genocols = ["name", "population", "locus", "genotype"]
        genocheck = intersect(genocols, names(loci))
        if genocheck != genocols
            throw(error("genodata missing columns $(symdiff(genocheck, genocols))"))
        end
        isempty(setdiff(meta.sampleinfo.name, loci.name.pool)) || throw(ArgumentError("metadata and genodata do not contain the same sample names"))
        isempty(setdiff(meta.sampleinfo.population, loci.population.pool)) || throw(ArgumentError("metadata and genodata do not contain the same population names"))
        new(meta, loci)
    end
end

PopData(data::DataFrame) = PopData(PopDataInfo(data), data)

# method to update PopDataInfo from PopData, all in one swoop
function PopDataInfo!(data::PopData)
    data.metadata.samples = length(data.genodata.name.pool)
    data.metadata.loci = length(data.genodata.locus.pool)
    data.metadata.populations = length(data.genodata.population.pool)
    filter!(:name => x -> x ∈ data.genodata.name.pool,  data.sampleinfo)
    filter!(:locus => x -> x ∈ data.genodata.locus.pool,  data.locusinfo)
    if "ploidy" ∈ names(data.sampleinfo)
        ploidy = unique(data.sampleinfo.ploidy)
        ploidy = length(ploidy) == 1 ? Int8(ploidy[1]) : Int8.(ploidy)
    else
        ploidy = Int8(0)
    end       
    data.metadata.ploidy = ploidy
    data.metadata.biallelic = data.metadata.biallelic ? true : isbiallelic(data)
    return
end

# method to update preexisting PopDataInfo with new genodata
# useful for getindex and creating new PopData from that
function PopDataInfo!(popdatainfo::PopDataInfo, genodata::DataFrame)
    popdatainfo.samples = length(genodata.name.pool)
    popdatainfo.loci = length(genodata.locus.pool)
    popdatainfo.populations = length(genodata.population.pool)
    filter!(:name => x -> x ∈ genodata.name.pool,  popdatainfo.sampleinfo)
    filter!(:locus => x -> x ∈ genodata.locus.pool,  popdatainfo.locusinfo)
    if "ploidy" ∈ names(popdatainfo.sampleinfo)
        ploidy = unique(popdatainfo.sampleinfo.ploidy)
        ploidy = length(ploidy) == 1 ? Int8(ploidy[1]) : Int8.(ploidy)
    else
        ploidy = Int8(0)
    end       
    popdatainfo.ploidy = ploidy
    popdatainfo.biallelic = popdatainfo.biallelic ? true : isbiallelic(genodata)
    return popdatainfo
end


"""
    Genotype::DataType
For convenience purposes, an alias for `NTuple{N, <:Signed} where N`, which is
the type describing individual genotypes in PopData. Specifically, there exist
`SNP` as an alias for `NTuple{N, Int8}` and `MSat` for `NTuple{N, Int16}`
"""
const Genotype = NTuple{N, <:Signed} where N

"""
    SNP::DataType
An alias for `NTuple{N, Int8}`
"""
const SNP = NTuple{N, Int8} where N
_SNP(geno) = isallmissing(geno) ? missing : SNP(geno)

"""
MSat::DataType
An alias for `NTuple{N, Int16}`
    """
const MSat = NTuple{N, Int16} where N
_MSat(geno) = isallmissing(geno) ? missing : MSat(geno)

"""
    GenoArray::DataType
An alias for an `AbstractVector` of elements `Missing`
and `Genotype`, which itself is of type `NTuple{N, <:Integer} where N`.
The definition as an `AbstractVector` adds flexibility for `SubArray`
cases.
"""
const GenoArray = AbstractVector{S} where S<:Union{Missing,Genotype}

function _ploidy2text(ploidy::Int8)
    ploidy == 0 ? "" :
    ploidy == 1 ? "Haploid" :
    ploidy == 2 ? "Diploid" :
    ploidy == 3 ? "Triploid" :
    ploidy == 4 ? "Tetraploid" :
    ploidy == 5 ? "Pentaploid" :
    ploidy == 6 ? "Hexaploid" : "Octaploid"
end

_ploidy2text(ploidy::Vector{Int8}) = "Mixed-ploidy"

function Base.show(io::IO, data::PopData)
    if occursin("Int16", string(eltype(data.genodata.genotype)))
        marker = "Microsatellite"
    else
        marker = "SNP"
    end
    ploidy = _ploidy2text(data.metadata.ploidy)
    ploidy = ploidy == "" ? "" : ploidy * ", "
    n_loc = data.metadata.loci
    println(io, "PopData", "{" * ploidy, n_loc, " " , marker * " loci}")
    println(io, "  Samples: $(data.metadata.samples)")
    print(io, "  Populations: $(data.metadata.populations)")
    allcols = vcat(names(data.sampleinfo), names(data.genodata)) |> unique
    extracols = symdiff(allcols, ["name", "population", "ploidy", "locus", "genotype"])
    if !isempty(extracols)
        print(io, "\n  Other Info: ", extracols)
    end
end

"""
    getindex(PopData, Symbol)
Return specific elements in a PopData objects, such at `sampleinfo`, `ploidy`, etc.
"""
function Base.getindex(data::PopData, idx::Symbol)
    if idx == :sampleinfo
        return data.sampleinfo
    elseif idx == :genodata
        return data.genodata
    elseif idx == :name
        return data.metadata..sampleinfo.name
    elseif idx == :population
        return data.sampleinfo.population
    elseif idx == :ploidy
        return data.sampleinfo.ploidy
    elseif idx == :locus
        return data.genodata.locus.pool
    elseif idx == :coordinates
        return data.sampleinfo[:, [:longitude, :latitude]]
    else
        throw(ArgumentError("Cannot directly index PopData with the \':$idx\' field"))
    end
end

"""
    getindex(PopData, args...)
    getindex(PopData, args..., cols)
Return a new PopData object or specific columns by indexing using standard DataFrames.jl conventions.
Expressions can be compounded using `.&` and `.|` operators.

**Example**
```julia
julia> cats = @nancycats ;
julia> cats[cats.genodata.name .== "N290"]
julia> cats[cats.genodata.name .∈ Ref(["N290", "N291]), :genotype]
```
"""
function Base.getindex(data::PopData, args)
    geno = getindex(data.genodata, args, :)
    transform!(
        geno,
        1 => (i -> PooledArray(i, compress = true)) => :name,
        2 => (i -> PooledArray(i, compress = true)) => :population,
        3 => (i -> PooledArray(i, compress = true)) => :locus,
        4
    )
    pdinfo = deepcopy(data.info)
    out = PopDataInfo!(pdinfo, geno)
    PopData(out, geno)
end
precompile(Base.getindex, (PopData, BitVector))


function Base.getindex(data::PopData, expression, cols)
    getindex(data.genodata, expression, cols)
end
precompile(Base.getindex, (PopData, BitVector, Colon))
precompile(Base.getindex, (PopData, BitVector, Symbol))
precompile(Base.getindex, (PopData, BitVector, Vector{Symbol}))
precompile(Base.getindex, (PopData, BitVector, Vector{Int64}))

"""
    getindex(PopData, NamedTuple)
    getindex(PopData, NamedTuple, cols)
Return a new PopData object or specific genodata columns by indexing a PopData object using a NamedTuple.
This method is syntactic sugar, but limited to basic `==` indexing.

**Example**
```julia
julia> cats = @nancycats ;
julia> a = (; locus = "fca8", population = ["1", "2", "3"]) ;
julia> cats[a]
PopData{Diploid, 1 Microsatellite loci}
  Samples: 44
  Populations: 3

julia> cats[a, :]
44×4 DataFrame
 Row │ name     population  locus   genotype   
     │ String7  String      String  Tuple…?    
─────┼─────────────────────────────────────────
   1 │ N215     1           fca8    missing
   2 │ N216     1           fca8    missing
   3 │ N217     1           fca8    (135, 143)
  ⋮  │    ⋮         ⋮         ⋮         ⋮
  42 │ N33      3           fca8    (123, 137)
  43 │ N34      3           fca8    (135, 139)
  44 │ N70      3           fca8    (143, 145)
                                38 rows omitted
```
"""
function Base.getindex(data::PopData, args::NamedTuple, cols)
    idx = mapreduce(.&, pairs(args)) do kv
        k = kv.first
        v = kv.second
        v isa AbstractVector ? data.genodata[:,k] .∈ Ref(v) : data.genodata[:,k] .== v
    end
    getindex(data.genodata, idx, cols)
end


function Base.getindex(data::PopData, args::NamedTuple)
    idx = mapreduce(.&, pairs(args)) do kv
        k = kv.first
        v = kv.second
        v isa AbstractVector ? data.genodata[:,k] .∈ Ref(v) : data.genodata[:,k] .== v
    end
    geno = getindex(data.genodata, idx, :)
    transform!(
        geno,
        1 => (i -> PooledArray(i, compress = true)) => :name,
        2 => (i -> PooledArray(i, compress = true)) => :population,
        3 => (i -> PooledArray(i, compress = true)) => :locus,
        4
    )
    pdinfo = deepcopy(data.info)
    out = PopDataInfo!(pdinfo, geno)
    PopData(out, geno)
end

function getindex_not(data::PopData, args::NamedTuple)
    idx = mapreduce(.&, pairs(args)) do kv
        k = kv.first
        v = kv.second
        v isa AbstractVector ? data.genodata[:,k] .∉ Ref(v) : data.genodata[:,k] .!= v
    end
    geno = getindex(data.genodata, idx, :)
    transform!(
        geno,
        1 => (i -> PooledArray(i, compress = true)) => :name,
        2 => (i -> PooledArray(i, compress = true)) => :population,
        3 => (i -> PooledArray(i, compress = true)) => :locus,
        4
    )
    pdinfo = deepcopy(data.info)
    out = PopDataInfo!(pdinfo, geno)
    PopData(out, geno)
end


# shortcut methods for convenience and less verbose typing
"""
    getproperty(data::PopData, field::Symbol)
    data.field

A convenience method to access certain elements in a `PopData` with fewer keystrokes. 
Essentially a standard `getproperty` call, except `sampleinfo` accesses `metadata.sampleinfo`,
`locusinfo` accesses `metadata.locusinfo`, and `info` is an alias for `metadata`.

## Example
```julia
cats = @nancycats ;
cats.metadata == cats.info
cats.metadata.sampleinfo == cats.sampleinfo
cats.metadata.locusinfo == cats.locusinfo
```
"""
function Base.getproperty(data::PopData, field::Symbol)
    if field == :sampleinfo
        return data.metadata.sampleinfo
    elseif field == :locusinfo
        return data.metadata.locusinfo
    elseif field == :info
        return data.metadata
    else
        return getfield(data, field)
    end
end