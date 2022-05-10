#=
General use utilities
=#

function Base.copy(data::PopData)
    PopData(deepcopy(data.metadata), copy(data.genodata))
end



function Base.size(data::PopData)
    return (samples = data.metadata.samples, loci = data.metadata.loci)
end

"""
    Base.sort(x::NTuple{N,T}) where N where T <: Signed 
Sort the integers within a Tuple and return the sorted Tuple.
"""
function Base.sort(x::NTuple{N,T}) where N where T <: Signed 
    Tuple(sort(SVector(x)))
end


"""
    convertcoord(coordinate::String)
Takes non-decimal-degree format as a `String` and returns it as a decimal degree
`Float32`. Can be broadcasted over an array of coordinate strings to convert them.
## Formatting requirements
- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as `missing` or the string `"missing"`
- Can mix colons and spaces (although it's bad practice)
### Example
```
julia> convertcoord("-41 31.52")
-41.5253f0
julia> convertcoord.(["-41 31.52", "25 11:54S"])
2-element Array{Float32,1}:
-41.5253
-25.1983
```
"""
function convertcoord(coordinate::String)
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

convertcoord(coordinate::Missing) = missing


"""
    dropmonomorphic(data::PopData; silent::Bool = false)
Return a `PopData` object omitting any monomorphic loci. Will inform you which
loci were removed.
"""
function dropmonomorphic(data::PopData; silent::Bool = false)
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samplenames(data)))
    monomorphs = string.(all_loci[[length(unique(i))==1 for i in prt]])
    if length(monomorphs) == 0
        #!silent && println("No monomorphic loci found.\n")
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
    dropmonomorphic!(data::PopData; silent::Bool = false)
Edit a `PopData` object in place by omitting any monomorphic loci. Will inform you which
loci were removed.
"""
function dropmonomorphic!(data::PopData; silent::Bool = false)
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samplenames(data)))
    monomorphs = string.(all_loci[[length(unique(i))==1 for i in prt]])
    if length(monomorphs) == 0
        #!silent && println("No monomorphic loci found.\n")
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
    dropmultiallelic(data::PopData)
Return a `PopData` object omitting loci that are not biallelic.
"""
function dropmultiallelic(data::PopData)
    if data.metadata.biallelic == true
        prinln("PopData object is already biallelic")
        return data
    end 
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samplenames(data)))
    nonbi = string.(all_loci[[!isbiallelic(i) for i in prt]])
    if length(nonbi) == 1
        @info "Removing 1 multiallelic locus"
        println()
    else
        @info "Removing $(length(nonbi)) multialleic loci"
        println()
    end
    _out = exclude(data, locus = nonbi)
    _out.metadata.biallelic = true
    return _out
end


"""
    dropmultiallelic!(data::PopData)
Edit a `PopData` object in place, removing loci that are not biallelic.
"""
function dropmultiallelic!(data::PopData)
    if data.metadata.biallelic == true
        println("PopData object is already biallelic\n")
        return data
    end 
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samplenames(data)))
    nonbi = string.(all_loci[[!isbiallelic(i) for i in prt]])
    if length(nonbi) == 1
        @info "Removing 1 multiallelic locus"
        println()
    else
        @info "Removing $(length(nonbi)) multialleic loci"
        println()
    end
    exclude!(data, locus = nonbi)
    data.metadata.biallelic = true
    return data
end

function truncatepath(text::String)
    width = displaysize(stdout)[2]
    if length(text) > width
        separated = split(text, "/")    
        newtxt = join(separated[[1,2,3]], "/") * "/…/" * join(separated[[end-1, end]],"/")
        return newtxt
    else
        return text   
    end
end
