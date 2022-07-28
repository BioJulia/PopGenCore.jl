"""
    allelecount(locus::T) where T<:GenoArray
Return the number of unique alleles present at a locus.
"""
@inline function allelecount(locus::AbstractVector{U})::Int64 where U<:Union{Missing, NTuple{N,T}} where N where T<: Union{Int8, Int16}
    out = 0
    uniq = T[]
    @inbounds for geno in skipmissing(locus)
        @inbounds @simd for allele in geno
            contains = allele ∉ uniq
            out += contains
            contains && push!(uniq, allele)
        end
    end
    return out
end

"""
    alleles(locus::T) where T<:GenoArray
Return an array of all the non-missing alleles of a locus.
"""
@inline function alleles(locus::AbstractVector{U}) where U<:Union{Missing, NTuple{N,T}} where N where T<: Union{Int8, Int16}
    skipm = skipmissing(locus)
    if isempty(skipm)
        return Vector{Union{Missing, T}}(undef, length(locus))
    else
        return T[j for i in skipm for j in i]
    end
end
precompile(alleles, (Vector{Union{Missing, NTuple{2, Int8}}},))
precompile(alleles, (Vector{Union{Missing, NTuple{2, Int16}}},))

"""
    alleles(locus::T, miss::Bool = false) where T<:GenoArray
Return an array of all the non-missing alleles of a locus. Use the second positional
argument as `true` to include missing values.
"""
@inline function alleles(locus::T, miss::Bool) where T<:GenoArray
    int_type = eltype(T) |> nonmissingtype |> eltype
    if isallmissing(locus)
        return Vector{Union{Missing, int_type}}(undef, length(locus))
    end
    alle_out = Union{Missing, int_type}[j for i in skipmissing(locus) for j in i]
    if miss == true
        nmiss = count(ismissing, locus)
        append!(alle_out, fill(missing, nmiss))
    end
    return alle_out
end
precompile(alleles, (Vector{Union{Missing, NTuple{2, Int8}}},Bool))
precompile(alleles, (Vector{Union{Missing, NTuple{2, Int16}}},Bool))


"""
    uniquealleles(locus::T) where T<:GenoArray
Return an array of all the unique non-missing alleles of a locus.
"""
@inline function uniquealleles(locus::AbstractVector{U})::Vector{T} where U<:Union{Missing, NTuple{N,T}} where N where T<: Union{Int8, Int16}
    out = T[]
    @inbounds for geno in skipmissing(locus)
        @inbounds @simd for allele in geno
            allele ∉ out && push!(out, allele)
        end
    end
    sort!(out)  # is this necessary?
    return out
end

"""
    uniquebialleles(locus::T) where T<:GenoArray
Return an array of all the unique non-missing alleles of a biallelic locus. This is similar
to `uniquealleles` but terminates once two alleles have been identified.
"""
@inline function uniquebialleles(locus::AbstractVector{U})::Vector{T} where U<:Union{Missing, NTuple{N,T}} where N where T<: Union{Int8, Int16}
    out = T[]
    @inbounds for geno in skipmissing(locus)
        @inbounds @simd for allele in geno
            allele ∉ out && push!(out, allele)
        end
        length(out) == 2 && break
    end
    sort!(out)  # is this necessary?
    return out
end

"""
    locidataframe(data::PopData)
Return a wide `DataFrame` of samples as columns, ommitting population information.

**Example**
```
julia> locidataframe(@nancycats)
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
function locidataframe(data::PopData)
    unstack(select(data.genodata, Not(:population)), :name, :genotype)[:, Not(:locus)]
end

#TODO make a SMatrix instead?
"""
    locimatrix(data::PopData)
Return a matrix of genotypes with dimensions `samples × loci`.
Rows are samples and columns are loci. Will return an error if ploidy varies between samples. 

**Example**

```julia
julia> locimatrix(@nancycats)
237×9 Matrix{Union{Missing, Tuple{Int16, Int16}}}:
 missing     …  (208, 208)
 missing        (208, 208)
 (135, 143)     (210, 210)
 (133, 135)     (208, 208)
 (133, 135)     (208, 208)
 (135, 143)  …  (208, 208)
 ⋮           ⋱
 (133, 141)     (208, 220)
 (133, 143)     (208, 208)
 (135, 141)     (208, 208)
 (137, 143)  …  (208, 208)
 (135, 141)     (208, 208)
```
"""
function locimatrix(data::PopData)
    dims = size(data)
    sort_df = issorted(data.genodata, [:name, :locus]) ? sort(data.genodata, [:name, :locus]) : data.genodata
    reshape(sort_df.genotype, (dims.samples, dims.loci)) |> collect
end


"""
    phasedmatrix(data::PopData)
Return a `Vector` of length `ploidy` composed of allele matrices with dimensions `samples × loci`.
Rows are samples and columns are loci. Will return an error if ploidy varies between samples. 

**Example**
```
julia> mtx = phasedmatrix(@nancycats)
2-element Array{Array{Union{Missing, Int16},2},1}:
 [missing 136 … 113 208; missing 146 … 113 208; … ; 137 130 … 113 208; 135 130 … missing 208]
 [missing 146 … 113 208; missing 146 … 113 208; … ; 143 136 … 117 208; 141 146 … missing 208]

julia> mtx[1]
237×9 Array{Union{Missing, Int16},2}:
    missing  136  139  116         156  142  199  113         208
    missing  146  139  120         156  142  185  113         208
 135         136  141  116         152  142  197  113         210
 133         138  139  116         150  142  199   91         208
 133         140  141  126         152  142  193  113         208
 135         136  145  120         150  148  193   91         208
 135         136  139  116         152  142  199  105         208
 135         136  135  120         154  142  193   91         208
 137         136  139  116         150  142  197  105         208
 135         132  141  120         150  148  197   91         208
 137         130  137  128         152  142  193   91         182
 129         130  135  126         144  140  193   91         182
   ⋮                                      ⋮                   
 133         136  135     missing  146  142  199  113         182
 133         136  135     missing  150  142  197  113         182
 133         130  141     missing  148  142  191     missing  208
 123         138  141     missing  148  142  191     missing  208
 123         138  139     missing  150  142  197     missing  208
 133         136  139     missing  150  142  197     missing  208
 133         130  139     missing  152  142  191     missing  208
 133         136  139     missing  150  142  199     missing  208
 133         130  135     missing  148  142  197     missing  208
 135         136  143     missing  144  142  191  113         208
 137         130  135     missing  150  142  193  113         208
 135         130  135     missing  150  142  197     missing  208
```
"""
function phasedmatrix(data::PopData)
    dims = size(data)
    ploidy = unique(data.sampleinfo.ploidy)
    ploidy = length(ploidy) != 1 ? error("Phasing will not work on mixed-ploidy samples") : ploidy[1]
    sort_df = issorted(data.genodata, [:name, :locus]) ? sort(data.genodata, [:name, :locus]) : data.genodata
    matrices = map(j -> map(i -> ismissing(i) ? missing : i[j] , sort_df.genotype), 1:ploidy)
    map(i -> collect(reshape(i, (dims.samples, dims.loci))), matrices)
end