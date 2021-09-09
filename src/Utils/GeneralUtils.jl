export copy, size, sort 
export generate_meta
export drop_monomorphic, drop_monomorphic!
export drop_multiallelic, drop_multiallelic!
export loci, samples, convert_coord

#=
General use utilities
=#


function Base.copy(data::PopData)
    PopData(copy(data.metadata), copy(data.genodata))
end

function Base.size(data::PopData)
    return (samples = size(data.metadata)[1], loci = length(loci(data)))
end

"""
    Base.sort(x::NTuple{N,T}) where N where T <: Signed 
Sort the integers within a Tuple and return the sorted Tuple.
"""
function Base.sort(x::NTuple{N,T}) where N where T <: Signed 
    Tuple(sort(SVector(x)))
end


function convert_coord(coordinate::String)
    lowercase(coordinate) == "missing" && return missing
    coord_strip = replace(uppercase(coordinate), r"[NSEW]" => "")
    split_coord = parse.(Float32, split(coord_strip, r"\s|:"))
    split_coord[2] /= 60.0
    if length(split_coord) == 3
        split_coord[3] /= 3600.0
    end
    conv = mapreduce(abs, +, split_coord)
    # N + E are positive | S + W are negative
    if split_coord[1] < 0 || occursin(r"[SW]", uppercase(coordinate))
        # negative
        return round(conv * -1, digits = 4)
    else
        # positive
        return round(conv, digits = 4)
    end
end


"""
    drop_monomorphic(data::PopData; silent::Bool = false)
Return a `PopData` object omitting any monomorphic loci. Will inform you which
loci were removed.
"""
function drop_monomorphic(data::PopData; silent::Bool = false)
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samples(data)))
    monomorphs = string.(all_loci[[length(unique(i))==1 for i in prt]])
    if length(monomorphs) == 0
        return data
    elseif !silent
        if length(monomorphs) == 1
            println("Removing monomorphic locus: " * monomorphs[1])
            println()
        else
            println("Removing $(length(monomorphs)) monomorphic loci:" * "\n $monomorphs")
            println()
        end
    end
    exclude(data, locus = loci_to_rm)
end


"""
    drop_monomorphic!(data::PopData; silent::Bool = false)
Edit a `PopData` object in place by omitting any monomorphic loci. Will inform you which
loci were removed.
"""
function drop_monomorphic!(data::PopData; silent::Bool = false)
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samples(data)))
    monomorphs = string.(all_loci[[length(unique(i))==1 for i in prt]])
    if length(monomorphs) == 0
        return data
    elseif !silent
        if length(monomorphs) == 1
            println("Removing monomorphic locus: " * monomorphs[1])
            println()
        else
            println("Removing $(length(monomorphs)) monomorphic loci:" * "\n $monomorphs")
            println()
        end
    end
    exclude!(data, locus = monomorphs)
end


"""
    drop_multiallelic(data::PopData)
Return a `PopData` object omitting loci that are not biallelic.
"""
function drop_multiallelic(data::PopData)
    if data.info.biallelic == true
        prinln("PopData object is already biallelic")
        return data
    end 
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samples(data)))
    nonbi = string.(all_loci[[!isbiallelic(i) for i in prt]])
    if length(nonbi) == 1
        @info "Removing 1 multiallelic locus"
        println()
    else
        @info "Removing $(length(nonbi)) multialleic loci"
        println()
    end
    exclude(data, locus = nonbi)
end


"""
    drop_multiallelic!(data::PopData)
Edit a `PopData` object in place, removing loci that are not biallelic.
"""
function drop_multiallelic!(data::PopData)
    if data.info.biallelic == true
        prinln("PopData object is already biallelic")
        return data
    end 
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samples(data)))
    nonbi = string.(all_loci[[!isbiallelic(i) for i in prt]])
    if length(nonbi) == 1
        @info "Removing 1 multiallelic locus"
        println()
    else
        @info "Removing $(length(nonbi)) multialleic loci"
        println()
    end
    exclude!(data, locus = nonbi)
    data.info.biallelic = true
    return data
end


"""
    generate_meta(data::DataFrame)
Given a genotype DataFrame formatted like `PopData.genodata`, generates a corresponding
`meta` DataFrame. In other words, it creates the `.metadata` part of `PopData` from the `.genodata` part.

**Example**
```
julia> cats = @nancycats ;

julia> cats_nometa = cats.genodata ;

julia> cats_meta = generate_meta(cats_nometa)
237×5 DataFrame
 Row │ name    population  ploidy  longitude  latitude 
     │ String  String      Int8    Float32?   Float32? 
─────┼─────────────────────────────────────────────────
   1 │ N215    1                2   missing   missing  
   2 │ N216    1                2   missing   missing  
   3 │ N217    1                2   missing   missing  
   4 │ N218    1                2   missing   missing  
   5 │ N219    1                2   missing   missing  
   6 │ N220    1                2   missing   missing  
   7 │ N221    1                2   missing   missing  
  ⋮  │   ⋮         ⋮         ⋮         ⋮         ⋮
 232 │ N295    17               2   missing   missing  
 233 │ N296    17               2   missing   missing  
 234 │ N297    17               2   missing   missing  
 235 │ N281    17               2   missing   missing  
 236 │ N289    17               2   missing   missing  
 237 │ N290    17               2   missing   missing  
                                       224 rows omitted
```
"""
function generate_meta(data::DataFrame)
    grp = groupby(data, :name)
    nms = map(z -> z.name, keys(grp))
    pops = [first(z.population) for z in grp]
    ploids = [find_ploidy(z.genotype) for z in grp]
    DataFrame(
        :name => nms,
        :population => pops,
        :ploidy => ploids,
        :longitude => Vector{Union{Missing, Float32}}(undef, (length(nms))),
        :latitude => Vector{Union{Missing, Float32}}(undef, (length(nms))),
        copycols = true
    )
end


"""
    loci(data::PopData)
Returns an array of strings of the loci names in a `PopData` object.
"""
function loci(data::PopData)
    unique(data.genodata.locus)
end

"""
    samples(data::PopData)
View individual/sample names in a `PopData`
"""
function samples(data::PopData)
    @view data.metadata[!, :name]
end

