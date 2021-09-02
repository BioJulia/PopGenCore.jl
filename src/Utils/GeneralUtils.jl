export copy, size, sort 
export generate_meta, loci_dataframe, loci_matrix
export drop_monomorphic, drop_monomorphic!
export drop_multiallelic, drop_multiallelic!
export loci, samples, convert_coord

#=
General use utilities
=#


function Base.copy(data::PopData)
    PopData(copy(data.metadatadata), copy(data.genodata))
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
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samples(data)))
    nonbi = string.(all_loci[[!isbiallelic(i) for i in prt]])
    if length(nonbi) == 0
        return data
    elseif length(nonbi) == 1
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
    all_loci = loci(data)
    prt = Base.Iterators.partition(data.genodata.genotype, length(samples(data)))
    nonbi = string.(all_loci[[!isbiallelic(i) for i in prt]])
    if length(nonbi) == 0
        return data
    elseif length(nonbi) == 1
        @info "Removing 1 multiallelic locus"
        println()
    else
        @info "Removing $(length(nonbi)) multialleic loci"
        println()
    end
    exclude!(data, locus = nonbi)
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
    loci_dataframe(data::PopData)
Return a wide `DataFrame` of samples as columns, ommitting population information.

**Example**
```
julia> loci_dataframe(@nancycats)
9×237 DataFrame. Omitted printing of 232 columns
│ Row │ N215       │ N216       │ N217       │ N218       │ N219       │
│     │ Tuple…?    │ Tuple…?    │ Tuple…?    │ Tuple…?    │ Tuple…?    │
├─────┼────────────┼────────────┼────────────┼────────────┼────────────┤
│ 1   │ missing    │ missing    │ (135, 143) │ (133, 135) │ (133, 135) │
│ 2   │ (136, 146) │ (146, 146) │ (136, 146) │ (138, 138) │ (140, 146) │
│ 3   │ (139, 139) │ (139, 145) │ (141, 141) │ (139, 141) │ (141, 145) │
│ 4   │ (116, 120) │ (120, 126) │ (116, 116) │ (116, 126) │ (126, 126) │
│ 5   │ (156, 156) │ (156, 156) │ (152, 156) │ (150, 150) │ (152, 152) │
│ 6   │ (142, 148) │ (142, 148) │ (142, 142) │ (142, 148) │ (142, 148) │
│ 7   │ (199, 199) │ (185, 199) │ (197, 197) │ (199, 199) │ (193, 199) │
│ 8   │ (113, 113) │ (113, 113) │ (113, 113) │ (91, 105)  │ (113, 113) │
│ 9   │ (208, 208) │ (208, 208) │ (210, 210) │ (208, 208) │ (208, 208) │
```
"""
function loci_dataframe(data::PopData)
    unstack(select(data.genodata, Not(:population)), :name, :genotype)[:, Not(:locus)]
end

#TODO make a SMatrix instead?
"""
    loci_matrix(data::PopData)
Return a matrix of genotypes with dimensions `samples × loci`.
Rows are samples and columns are loci. Will return an error if ploidy varies between samples. 

**Example**
```
julia> loci_matrix(@nancycats)
237×9 Array{Union{Missing, Tuple{Int16,Int16}},2}:
 missing     (136, 146)  (139, 139)  …  (199, 199)  (113, 113)  (208, 208)
 missing     (146, 146)  (139, 145)     (185, 199)  (113, 113)  (208, 208)
 (135, 143)  (136, 146)  (141, 141)     (197, 197)  (113, 113)  (210, 210)
 (133, 135)  (138, 138)  (139, 141)     (199, 199)  (91, 105)   (208, 208)
 (133, 135)  (140, 146)  (141, 145)     (193, 199)  (113, 113)  (208, 208)
 (135, 143)  (136, 146)  (145, 149)  …  (193, 195)  (91, 113)   (208, 208)
 (135, 135)  (136, 146)  (139, 145)     (199, 199)  (105, 113)  (208, 208)
 (135, 143)  (136, 146)  (135, 149)     (193, 197)  (91, 91)    (208, 212)
 (137, 143)  (136, 146)  (139, 139)     (197, 197)  (105, 113)  (208, 212)
 (135, 135)  (132, 132)  (141, 145)     (197, 197)  (91, 105)   (208, 208)
 (137, 141)  (130, 136)  (137, 145)  …  (193, 199)  (91, 91)    (182, 182)
 (129, 133)  (130, 136)  (135, 145)     (193, 199)  (91, 113)   (182, 208)
 ⋮                                   ⋱                          
 (133, 135)  (136, 136)  (135, 139)  …  (199, 199)  (113, 113)  (182, 182)
 (133, 141)  (136, 136)  (135, 139)     (197, 197)  (113, 113)  (182, 208)
 (133, 141)  (130, 146)  (141, 141)     (191, 199)  missing     (208, 208)
 (123, 133)  (138, 138)  (141, 145)     (191, 197)  missing     (208, 208)
 (123, 133)  (138, 138)  (139, 139)     (197, 199)  missing     (208, 208)
 (133, 141)  (136, 146)  (139, 139)  …  (197, 197)  missing     (208, 208)
 (133, 141)  (130, 136)  (139, 145)     (191, 199)  missing     (208, 208)
 (133, 141)  (136, 146)  (139, 145)     (199, 199)  missing     (208, 220)
 (133, 143)  (130, 130)  (135, 145)     (197, 197)  missing     (208, 208)
 (135, 141)  (136, 144)  (143, 143)     (191, 197)  (113, 117)  (208, 208)
 (137, 143)  (130, 136)  (135, 145)  …  (193, 199)  (113, 117)  (208, 208)
 (135, 141)  (130, 146)  (135, 139)     (197, 197)  missing     (208, 208)
 ```
"""
function loci_matrix(data::PopData)
    dims = size(data)
    sort_df = issorted(data.genodata, [:name, :locus]) ? sort(data.genodata, [:name, :locus]) : data.genodata
    reshape(sort_df.genotype, (dims.samples, dims.loci)) |> collect
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

