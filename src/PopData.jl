export PopObj, PopData, show, Genotype, GenoArray, SNP, MSat
export metadata, genodata, getindex

"""
    AbstractType PopObj
Generic AbstractType for use in PopGen.jl
"""
abstract type PopObj end

"""
```
PopData
    metadata::DataFrame
    genodata::DataFrame
```
The data struct used for the PopGen population genetics ecosystem. You are
**strongly** discouraged from manually creating tables to pass into a PopData,
and instead should use the provided file importers and utilities.

- `metadata` DataFrame of sample information with the columns:
    - `name` - sample names [`String`]
    - `population` - population names [`String`]
    - `ploidy` - sample ploidy [`Int8`]
    - `longitude` - longitude values [`Float32`, optional]
    - `latitude` - latitude values   [`Float32`, optional]
- `genodata` DataFrame of sample genotype records
    - `name` - the individual/sample names [`PooledArray`]
    - `population` - population names [`PooledArray`]
    - `locus` - locus names [`PooledArray`]
    - `genotype` - genotype values [`NTuple{N,Signed}`]
"""
struct PopData <: PopObj
    metadata::DataFrame
    genodata::DataFrame
    function PopData(meta::DataFrame, loci::DataFrame)
        metacols = ["name", "population", "ploidy"]
        metacheck = intersect(metacols, names(meta))
        if metacheck != metacols
            throw(error("metadata missing columns $(symdiff(metacheck, metacols))"))
        end
        genocols = ["name", "population", "locus", "genotype"]
        genocheck = intersect(genocols, names(loci))
        if genocheck != genocols
            throw(error("genodata missing columns $(symdiff(genocheck, genocols))"))
        end
        if !issorted(loci, [:locus, :population, :name], lt = natural)
            sort!(loci, [:locus, :population, :name], lt = natural)
        end
        sort(meta.name) != sort(loci.name.pool) && throw(ArgumentError("meta and loci dataframes do not contain the same sample names"))
        sort(unique(meta.population)) != sort(loci.population.pool) && throw(ArgumentError("metadata and genotypes dataframes do not contain the same population names"))
        new(meta, loci)
    end
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


function Base.show(io::IO, data::PopData)
    if occursin("Int16", string(eltype(data.genodata.genotype)))
        marker = "Microsatellite"
    else
        marker = "SNP"
    end
    if "ploidy" ∈ names(data.metadata)
        ploidy = unique(data.metadata.ploidy) |> sort
        if length(ploidy) == 1
            ploidy = first(ploidy)
            ploidytext = 
                ploidy == 1 ? "Haploid" :
                ploidy == 2 ? "Diploid" :
                ploidy == 3 ? "Triploid" :
                ploidy == 4 ? "Tetraploid" :
                ploidy == 5 ? "Pentaploid" :
                ploidy == 6 ? "Hexaploid" : "Octaploid"
        else
            ploidytext = "Mixed-ploidy"
        end
    else
        ploidytext = "Unknown-ploidy"
    end
    n_loc = length(loci(data))
    println(io, "PopData", "{" * ploidytext * ", ", n_loc, " " , marker * " loci}")
    println(io, "  Samples: $(length(data.metadata.name))")
    print(io, "  Populations: $(length(data.genodata.population.pool))")
    if "longitude" ∈ names(data.metadata)
        miss_count = count(ismissing, data.metadata.longitude)
        if miss_count == length(data.metadata.longitude)
            print("")
        elseif iszero(miss_count)
            print(io, "\n  Coordinates: present")
        else
            print(io, "\n  Coordinates: present (", count(ismissing, data.metadata.longitude), " missing)")
        end
    end
    allcols = vcat(names(data.metadata), names(data.genodata)) |> unique
    extracols = symdiff(allcols, ["name", "population", "ploidy", "longitude", "latitude", "locus", "genotype"])
    if !isempty(extracols)
        print(io, "\n  Other Info: ", extracols)
    end
end


"""
    metadata(::PopData)
Method to show the `PopData` `metadata` field. 
"""
function metadata(data::PopData)
    s,f  = size(data.metadata)
    dimtext = "(" * string(s) * " samples, " * string(f) * " fields)"
    show(data.metadata, show_row_number = false, title = "Metadata of PopData " * dimtext )
end

"""
    genodata(::PopData)
Method to show the `PopData` `genodata` field. 
"""
function genodata(data::PopData) 
    l = length(data.genodata.locus.pool)
    s = length(data.metadata.name)
    dimtext = "(" * string(s) * " samples, " * string(l) * " loci)"
    show(data.genodata, show_row_number = false, title = "Genotype information of PopData " * dimtext )
end


function Base.getindex(data::PopData, idx::Symbol)
    if idx == :metadata
        return data.metadata
    elseif idx == :genodata
        return data.genodata
    elseif idx == :name
        return data.metadata.name
    elseif idx == :population
        return data.metadata.population
    elseif idx == :ploidy
        return data.metadata.ploidy
    elseif idx == :locus
        return data.genodata.locus.pool
    elseif idx == :coordinates
        return data.metadata[:, [:longitude, :latitude]]
    else
        throw(ArgumentError("Cannot index PopData with the \':$idx\' field"))
    end
end

function Base.getindex(data::PopData, args...)
    getindex(data.genodata, args...)
end
