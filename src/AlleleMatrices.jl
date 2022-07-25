"""
    matrix(data::PopData, matrixtype::Union{String, Symbol} = "count", missings = "mean", scale = false, center = false)
Return a matrix of allele counts or frequencies per genotype where rows are samples
and columns are the occurence count or frequency of an allele for that locus in that sample.
Loci and alleles are sorted alphanumerically. Setting `scale` or `center` as `true` will
compute allele frequencies regardless of the `by` keyword.

### Positional Arguments
- `data`: a PopData object
- `matrixtype`: a `String` or `Symbol` of `count` or `frequency` (default: `frequency`)

### Keyword Arguments
- `missings`: a `String` denoting how to handle missing values when outputting `frequency` (default: `mean`)
    - `"missing"`: fallback method to keep `missing` values as they are
    - `"zero"`: replace `missing` values with `0`
    - `"mean"`: replace `missing` values with the mean frequency for that allele in that locus
- `scale`: a 'Bool' of whether to z-score scale allele frequencies (default: `false`)
- `center`: a `Bool` of whether to center the allele frequencies (default: 'false')

**Example**

```
julia> cats = @nancycats ;

julia> cnts = matrix(cats, "count") ;  cnts[1:5,1:6]
5×6 Matrix{Union{Missing, Int8}}:
  missing   missing   missing   missing   missing   missing
  missing   missing   missing   missing   missing   missing
 0         0         0         0         0         0
 0         0         0         0         0         0
 0         0         0         0         0         0

julia> frq = matrix(cats, "frequency") ;  frq[1:5,1:6]
5×6 Matrix{Union{Missing, Float32}}:
 0.00460829  0.00460829  0.0276498  0.133641  0.00460829  0.0921659
 0.00460829  0.00460829  0.0276498  0.133641  0.00460829  0.0921659
 0.0         0.0         0.0        0.0       0.0         0.0
 0.0         0.0         0.0        0.0       0.0         0.0
 0.0         0.0         0.0        0.0       0.0         0.0

 julia> frq = matrix(cats, "frequency", missings = "zero") ;  frq[1:5,1:6]
 5×6 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0

 julia> frq = matrix(cats, missings = "mean", scale = true, center = true) ;  frq[1:5,1:6]
 5×6 Matrix{Float32}:
  7.17017f-9   7.17017f-9   0.0        0.0        7.17017f-9   0.0
  7.17017f-9   7.17017f-9   0.0        0.0        7.17017f-9   0.0
 -0.0709577   -0.0709577   -0.175857  -0.394198  -0.0709577   -0.300797
 -0.0709577   -0.0709577   -0.175857  -0.394198  -0.0709577   -0.300797
 -0.0709577   -0.0709577   -0.175857  -0.394198  -0.0709577   -0.300797
```
 """
function matrix(data::PopData, matrixtype::Union{String, Symbol} = "frequency"; missings::String = "mean", scale::Bool = false, center::Bool = false)
    mthd = string(matrixtype)
    mthd ∉  ["frequency", "count"] && throw(ArgumentError("Matrix type $matrixtype not recognized. Choose from either \"count\" or \"frequency\" methods"))
    mthd = matrixtype == "frequency" ? matrixtype * missings : matrixtype
    mtx = _matrix(data, Val(Symbol(mthd)))    
    if (mthd != "count") & (scale | center)
        mtx = standardize(ZScoreTransform, mtx, dims = 1, scale = scale, center = center)
        # replace almost-zero values caused by missing values with 0.0
        replace!(x ->  0 < x < (10^-9) ? 0.0 : x, mtx)
    end
    return mtx
end

function _setcounts(q, r)
    l = 0
    @inbounds for i in r
        l += length(i)
    end
    cnt = Vector{eltype(eltype(r))}(undef, l)
    idx = 0
    @inbounds for (i,j) in enumerate(r)
        @inbounds geno = q[i]
        @inbounds @simd for h in j 
            idx += 1
            @inbounds cnt[idx] = geno === missing ? -1 : count(==(h), geno)
        end
    end
    return cnt
end


"""
    _matrix(data::PopData, ::Val{:count})
Create a matrix of allele count per genotype where rows are samples
and columns are the occurence count of an allele for that locus in that sample.
Missing values are preserved as `-1`.
"""
function _matrix(data::PopData, ::Val{:count})
    gmtx = locimatrix(data)
    allalleles = Tuple(uniquealleles(i) for i in eachcol(gmtx))
    mapreduce(hcat, eachrow(gmtx)) do smple
        _setcounts(smple, allalleles)
    end |> permutedims
end
precompile(_matrix, (PopData,Val{:count}))

"""
    _matrix(data::PopData, ::Val{:frequencyzero})
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are replaced by zeros.
"""
function _matrix(data::PopData, ::Val{:frequencyzero})
    out = matrix(data, Val(:count))
    replace!(out, -1 => 0)
    out ./ data.sampleinfo.ploidy
end
precompile(_matrix, (PopData,Val{:frequencyzero}))


"""
    _matrix(data::PopData, ::Val{:frequencymean})
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are replaced by the global mean allele frequency.
"""
function _matrix(data::PopData, ::Val{:frequencymean})
    counts = @inbounds matrix(data, Val(:count)) ./ data.sampleinfo.ploidy
    map(eachcol(counts)) do alcol
        colmean = mean([x for x in alcol if x >= 0])
        replace!(x -> x < 0 ? colmean : x, alcol)
    end
    return counts
end
precompile(_matrix, (PopData,Val{:frequencymean}))


"""
    _matrix(data::PopData, ::Val{:frequencymissing})
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are kept as `missing`.
"""
function _matrix(data::PopData, ::Val{:frequencymissing})
    out = allowmissing(matrix(data, Val(:count)))
    replace!(out, -1 => missing)
    out ./ data.sampleinfo.ploidy
end
precompile(_matrix, (PopData,Val{:frequencymissing}))


"""
    featurematrix(data::PopData, matrixtype::Union{String, Symbol} = "genotype")

### Positional Arguments
    - `data`: a PopData object
    - `matrixtype`: a `String` or `Symbol` of `genotype`, or `allele` (default: `genotype`)

**genotype feature matrix**

Return a matrix of dummy-encoded genotypes (0,1,2...), where rows correspond with samples and columns correspond to loci.
Missing genotypes are encoded as `-1`. For biallelic loci, `0` encodes homozygous for allele 1, `1` encodes for a heterozygote,
and `2` encodes for homozygous allele 2.

**allele feature matrix**

Return a matrix of dummy-encoded alleles (0,1), where rows correspond with samples and columns correspond to alleles within loci, such
that there are as many columns per locus as alleles for that locus. Missing alleles (from missing genotypes) are encoded as `-1`.

**Example**
```
julia> cats = @nancycats ;
julia> featurematrix(cats)
237×9 Matrix{Int8}:
 -1   0   0   0   0  0   0   0   0
 -1   1   1   1   0  0   1   0   0
  0   0   2   2   1  1   2   0   1
  1   2   3   3   2  0   0   1   0
  1   3   4   4   3  0   3   0   0
  ⋮                  ⋮
 49   0   1  -1  36  0   0  -1  13
 48   6   8  -1  25  1   2  -1   0
 29   9  23  -1   7  3  26  14   0
  3   5   8  -1   2  1   3  14   0
 29  10  16  -1   4  3   2  -1   0

 julia> featurematrix(cats, "allele")
237×108 Matrix{Int8}:
 -1  -1  -1  -1  -1  …  0  0  0  0  0  0
 -1  -1  -1  -1  -1     0  0  0  0  0  0
  0   0   0   0   0     0  0  0  0  0  0
  0   0   0   0   0     0  0  0  0  0  0
  0   0   0   0   0     0  0  0  0  0  0
  ⋮                  ⋱           ⋮
  0   0   0   0   0     0  0  0  1  0  0
  0   0   0   0   0     0  0  0  0  0  0
  0   0   0   0   0     0  0  0  0  0  0
  0   0   0   0   0  …  0  0  0  0  0  0
  0   0   0   0   0     0  0  0  0  0  0 
```
"""
function featurematrix(data::PopData, matrixtype::Union{String, Symbol} = "genotype")::Matrix{Int8}
    string(matrixtype) ∉ ["genotype", "allele"] && throw(ArgumentError("Matrix type $matrixtype not recognized. Please use one of \"genotype\" or \"allele\"."))
    _featurematrix(data, Val(Symbol(matrixtype))) 
end

"""
    _featurematrix(data::PopData, ::Val{:genotype})
Return a matrix of dummy-encoded genotypes (0,1,2...), where rows correspond with samples and columns correspond to loci.
Missing genotypes are encoded as `-1`. For biallelic loci, `0` encodes homozygous for allele 1, `1` encodes for a heterozygote,
and `2` encodes for homozygous allele 2.
"""
function _featurematrix(data::PopData, ::Val{:genotype})::Matrix{Int8}
    genomtx = locimatrix(data)
    if isbiallelic(data)
        mapreduce(hcat, eachcol(genomtx)) do i
            uniq_genos = unique(skipmissing(i))
            # find the heterozygous genotype
            het = findfirst(ishet, uniq_genos)
            # get the indices for the other ones
            others = setdiff(eachindex(uniq_genos), het)
            # setup the new dummy variables as 0:n_alleles-1 excecpt omitting 1
            newvals = Int8.(0:length(others))
            newvals = newvals[newvals .!= 1]
            gdict = Dict(zip(uniq_genos[others], newvals))
            # add the heterozygote at the end
            gdict[uniq_genos[het]] = Int8(1)
            #gdict = Dict(zip(uniq_genos, Int8.(0:length(uniq_genos)-1)))
            Int8[get(gdict, geno, Int8(-1)) for geno in i]
        end
    else
        mapreduce(hcat, eachcol(genomtx)) do i
            uniq_genos = unique(skipmissing(i))
            gdict = Dict(zip(uniq_genos, Int8.(0:length(uniq_genos)-1)))
            Int8[get(gdict, geno, Int8(-1)) for geno in i]
        end
    end
end
precompile(_featurematrix, (PopData,Val{:genotype}))

"""
    _featurematrix(data::PopData, ::Val{:allele})
Return a matrix of dummy-encoded alleles (0,1,2...), where rows correspond with samples and columns correspond to loci.
Missing genotypes are encoded as `-1`.
"""
function _featurematrix(data::PopData, ::Val{:allele})::Matrix{Int8}
    allecounts = matrix(data, Val(:count))
    replace!(allecounts) do i
        i > 0 ? Int8(1) : i
    end
    return allecounts
end
precompile(_featurematrix, (PopData, Val{:allele}))