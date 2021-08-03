using DataFrames, PooledArrays

struct PopData2
    metadata::DataFrame
    genotypes::DataFrame
    samples::Vector{String}
    populations::Vector{String}
    loci::Vector{String}
end

x = @nancycats

pd2pd2(x::PopData) = PopData2(x.meta, x.loci, x.loci.name.pool, x.loci.population.pool, x.loci.locus.pool)

function Base.show(io::IO, data::PopData2)
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
    n_loc = length(data.loci)
    println(io, "PopData", "{" * ploidytext * ", ", n_loc, " " , marker * " loci}")
    print(io, "  Samples: ") ; printstyled(io, length(data.samples), "\n", bold = true)
    print(io, "  Populations: ") ; printstyled(io, length(data.populations), bold = true)
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

## Notes

1. :name changing is bidirectional. Meaning you 
can modify the name of a sample in .metadata and it will replace all occurences of it in .genotypes

2. you can modify PopData.samples, PopData.loci, PopData.populations and it will apply all those
changes to both dataframes

3. adding a new population in .metadata will break the sync with the pools


