export PopObj, PopData, show, Genotype, GenoArray

"""
    AbstractType PopObj
Generic AbstractType for use in PopGen.jl
"""
abstract type PopObj end

"""
```
PopData
    metadata::DataFrame
    genotypes::DataFrame
```
The data struct used for the PopGen population genetics ecosystem. You are
**strongly** discouraged from manually creating tables to pass into a PopData,
and instead should use the provided file importers.

- `metadata` ::DataFrame individual/sample data with the columns:
    - `name` ::String the individual/sample names
    - `population` ::String population names
    - `ploidy` ::Int8 ploidy in order of `ind`
    - `longitude` ::Float64 longitude values
    - `latitude` ::Float64 latitude values
- `genotypes` ::DataFrame Long-format table of sample genotype records
    - name` ::PooledArray the individual/sample names
    - `population`::PooledArray of population names
    - `locus` ::PooledArray of locus names
    - `genotype` Tuple of Int8 or Int16 depending on SNP or microsatellite
"""
struct PopData <: PopObj
    metadata::DataFrame
    genotypes::DataFrame
    function PopData(meta::DataFrame, loci::DataFrame)
        sort!(loci, [:locus, :population, :name], lt = natural)
        sort(meta.name) != sort(loci.name.pool) && throw(ArgumentError("meta and loci dataframes do not contain the same sample names"))
        sort(unique(meta.population)) != sort(loci.population.pool) && throw(ArgumentError("metadata and genotypes dataframes do not contain the same population names"))
        new(meta, loci)
    end
end


"""
    Genotype::DataType
For convenience purposes, an alias for `NTuple{N, <:Integer} where N`, which is
the type describing individual genotypes in PopData. Specifically, there exist
`SNP` as an alias for `NTuple{N, Int8}` and `MSat` for `NTuple{N, Int16}`
"""
const Genotype = NTuple{N, <:Integer} where N

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
For convenience purposes, an alias for an `AbstractVector` of elements `Missing`
and `Genotype`, which itself is of type `NTuple{N, <:Integer} where N`.
The definition as an `AbstractVector` adds flexibility for `SubArray`
cases.
"""
const GenoArray = AbstractVector{S} where S<:Union{Missing,Genotype}


function Base.show(io::IO, data::PopData)
    if occursin("Int16", string(eltype(data.genotypes.genotype)))
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
    n_loc = length(data.genotypes.locus.pool)
    println(io, "PopData", "{" * ploidytext * ", ", n_loc, " " , marker * " loci}")
    println(io, "  Samples: $(length(data.metadata.name))") #; printstyled(io, length(data.samples), "\n", bold = true)
    print(io, "  Populations: $(length(data.genotypes.population.pool))") # ; printstyled(io, length(data.populations), bold = true)
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
    allcols = vcat(names(data.metadata), names(data.genotypes)) |> unique
    extracols = symdiff(allcols, ["name", "population", "ploidy", "longitude", "latitude", "locus", "genotype"])
    if !isempty(extracols)
        print(io, "\n  Other Info: ", extracols)
    end
end

#=
function metadata(data::PopData) 
    s,f  = size(data.metadata)
    dimtext = "(" * string(s) * " samples, " * string(f) * " fields)"
    show(data.metadata, show_row_number = false, title = "Metadata of PopData " * dimtext )
end

function genotypes(data::PopData) 
    l = length(data.genotypes.locus.pool)
    s = length(data.metadata.name)
    dimtext = "(" * string(s) * " samples, " * string(l) * " loci)"
    show(data.genotypes, show_row_number = false, title = "Genotype information of PopData " * dimtext )
end
=#