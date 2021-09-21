export PopObj, PopData, PopDataInfo, show, Genotype, GenoArray, SNP, MSat
export metadata, genodata, getindex

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
The data struct used internally as `PopData.info` fields to store basic information
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
    sampleinfo = unique(dropmissing(genodf), :name)
    insertcols!(sampleinfo, :ploidy => length.(sampleinfo.genotype))
    select!(sampleinfo, :name, :population, :ploidy)
    ploidy = unique(sampleinfo.ploidy)
    ploidy = length(ploidy) == 1 ? Int8(ploidy[1]) : Int8.(ploidy)
    PopDataInfo(
        length(genodf.name.pool),
        sampleinfo,
        length(genodf.locus.pool),
        DataFrame(:chromosome => Int8(0), :locus => genodf.locus.pool, :cm => Int8(0), :bp => Int8(0)),
        length(genodf.population.pool),
        ploidy,
        isbiallelic(genodf)
    )
end

function Base.show(io::IO, data::PopDataInfo)
    println(io, " ploidy:        ", data.ploidy...)
    println(io, " # loci:        ", data.loci)
    println(io, " # samples:     ", data.samples)
    println(io, " # populations: ", data.populations)
    println(io, " biallelic:     ", data.biallelic)
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
        sort(meta.sampleinfo.name) != sort(loci.name.pool) && throw(ArgumentError("metadata and genodata do not contain the same sample names"))
        sort(unique(meta.sampleinfo.population)) != sort(loci.population.pool) && throw(ArgumentError("metadata and genodata do not contain the same population names"))
        new(PopDataInfo(loci), loci)
    end
end

PopData(data::DataFrame) = PopData(PopDataInfo(data), data)

function PopDataInfo!(data::PopData)
    data.metadata.samples = length(data.genodata.name.pool)
    data.metadata.loci = length(data.genodata.locus.pool)
    data.metadata.populations = length(data.genodata.population.pool)
    filter!(:name => x -> x .== data.genodata.name.pool,  data.metadata.sampleinfo)
    filter!(:locus => x -> x .== data.genodata.locus.pool,  data.metadata.locusinfo)
    if "ploidy" ∈ names(data.metadata.sampleinfo)
        ploidy = unique(data.metadata.sampleinfo.ploidy)
        ploidy = length(ploidy) == 1 ? ploidy[1] : ploidy
    else
        ploidy = Int8(0)
    end       
    data.metadata.ploidy = ploidy
    data.metadata.biallelic = data.metadata.biallelic ? true : isbiallelic(data)
    return
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
_SNP(geno) = all(ismissing.(geno)) ? missing : SNP(geno)

"""
MSat::DataType
An alias for `NTuple{N, Int16}`
    """
const MSat = NTuple{N, Int16} where N
_MSat(geno) = all(ismissing.(geno)) ? missing : MSat(geno)

"""
    GenoArray::DataType
An alias for an `AbstractVector` of elements `Missing`
and `Genotype`, which itself is of type `NTuple{N, <:Integer} where N`.
The definition as an `AbstractVector` adds flexibility for `SubArray`
cases.
"""
const GenoArray = AbstractVector{S} where S<:Union{Missing,Genotype}
#const SnpArray = PooledVector{Union{Missing, Tuple{Int8, Int8}}, UInt8, Vector{UInt8}} 

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
    if "longitude" ∈ names(data.metadata.sampleinfo)
        miss_count = count(ismissing, data.metadata.sampleinfo.longitude)
        if miss_count == length(data.metadata.sampleinfo.longitude)
            print("")
        else
            print(io, "\n  Coordinates: present")
        end
    end
    allcols = vcat(names(data.metadata.sampleinfo), names(data.genodata)) |> unique
    extracols = symdiff(allcols, ["name", "population", "ploidy", "locus", "genotype"])
    if !isempty(extracols)
        print(io, "\n  Other Info: ", extracols)
    end
end


"""
    sampleinfo(::PopData)
Method to show the `PopData` `metadata` field. 
"""
function sampleinfo(data::PopData)
    s,f  = size(data.metadata.sampleinfo)
    dimtext = "(" * string(s) * " samples, " * string(f) * " fields)"
    show(data.metadata.sampleinfo, show_row_number = false, title = "Sample Information of PopData " * dimtext )
end

"""
    genodata(::PopData)
Method to show the `PopData` `genodata` field. 
"""
function genodata(data::PopData) 
    l = length(data.genodata.locus.pool)
    s = length(data.metadata.sampleinfo.name)
    dimtext = "(" * string(s) * " samples, " * string(l) * " loci)"
    show(data.genodata, show_row_number = false, title = "Genotype information of PopData " * dimtext )
end

function locusinfo(data::PopData)
    s,f  = size(data.metadata.locusinfo)
    dimtext = "(" * string(s) * " loci, " * string(f) * " fields)"
    show(data.metadata.locusinfo, show_row_number = false, title = "Locus Information of PopData " * dimtext )
end

function Base.getindex(data::PopData, idx::Symbol)
    if idx == :sampleinfo
        return data.metadata.sampleinfo
    elseif idx == :genodata
        return data.genodata
    elseif idx == :name
        return data.metadata..sampleinfo.name
    elseif idx == :population
        return data.metadata.sampleinfo.population
    elseif idx == :ploidy
        return data.metadata.sampleinfo.ploidy
    elseif idx == :locus
        return data.genodata.locus.pool
    elseif idx == :coordinates
        return data.metadata.sampleinfo[:, [:longitude, :latitude]]
    else
        throw(ArgumentError("Cannot directly index PopData with the \':$idx\' field"))
    end
end

function Base.getindex(data::PopData, args...)
    geno = getindex(data.genodata, args...)
    transform!(
        geno,
        1 => (i -> PooledArray(i, compress = true)) => :name,
        2 => (i -> PooledArray(i, compress = true)) => :population,
        3 => (i -> PooledArray(i, compress = true)) => :locus,
        4
    )
    #newmeta = intersect(data.metadata.sampleinfo.name, geno.name.pool) != data.metadata.sampleinfo.name ? data.metadata.sampleinfo[data.metadata.sampleinfo.name .∈ Ref(geno.name.pool), :] : data.metadata.sampleinfo
    PopData(geno)
end