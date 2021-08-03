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
    println(io, "PopData Object")
    if occursin("Int16", string(eltype(data.genotypes.genotype)))
        marker = "Microsatellite"
    else
        marker = "SNP"
    end
    n_loc = length(data.loci)
    printstyled(io, n_loc, " ", marker, " loci", "\n" , bold = true)
    if "ploidy" ∈ names(data.metadata)
        ploidy = unique(data.metadata.ploidy) |> sort
        if length(ploidy) == 1
            print(io, "  Ploidy: ") ; printstyled(io, ploidy |> join, "\n", bold = true)
        else
            print(io, "  Ploidy (varies): ")
            printstyled(io, ploidy[1], bold = true); [printstyled(io, ", $i", bold = true) for i in ploidy[2:end]]
            print(io, "\n")
        end
    else
        print(io, "  Ploidy:") ; printstyled(io, " absent\n", color = :yellow)
    end
    print(io, "  Samples: ") ; printstyled(io, length(data.samples), "\n", bold = true)
    print(io, "  Populations: ") ; printstyled(io, length(data.populations), bold = true)

    if "longitude" ∈ names(data.metadata)
        miss_count = count(ismissing, data.metadata.longitude)
        if miss_count == length(data.metadata.longitude)
            return
        elseif iszero(miss_count)
            print(io, "\n  Coordinates: present")
        else
            print(io, "\n  Coordinates: present (", count(ismissing, data.metadata.longitude), " missing)")
        end
    end
    # add extracols section
end
