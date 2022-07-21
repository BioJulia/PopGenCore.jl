"""
    matrix(data::PopData, matrixtype::String = "count", missings = "mean", scale = false, center = false)
Return a matrix of allele counts or frequencies per genotype where rows are samples
and columns are the occurence count or frequency of an allele for that locus in that sample.
Loci and alleles are sorted alphanumerically. Setting `scale` or `center` as `true` will
compute allele frequencies regardless of the `by` keyword.

### Positional Arguments
- `data`: a PopData object
- `matrixtype`: a `String` of `count` or `frequency` (default: `frequency`)

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
function matrix(data::PopData, matrixtype::String = "frequency"; missings::String = "mean", scale::Bool = false, center::Bool = false)
    if occursin(lowercase(matrixtype), "counts") && !any([scale, center])
        return countmatrix(data)
    elseif occursin(lowercase(matrixtype), "frequency") 
        freqs = 
            if lowercase(missings) == "mean"
                freqmatrix_mean(data)
            elseif lowercase(missings) == "zero"
                freqmatrix_zero(data)
            elseif lowercase(missings) == "missing"
                freqmatrix_missing(data)
            else
                throw(ArgumentError("use one of \"zero\" or \"mean\" for handling missing values"))
            end
        if any([scale, center])
            return freqmatrix_scale(freqs, scale, center)
        else
            return freqs
        end
    else
        throw(ArgumentError("Choose from either \"count\" or \"frequency\" methods"))
    end
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
    countmatrix(data::PopData)
Create a matrix of allele count per genotype where rows are samples
and columns are the occurence count of an allele for that locus in that sample.
Missing values are preserved as `-1`.
"""
function countmatrix(data::PopData)
    gmtx = locimatrix(data)
    allalleles = Tuple(uniquealleles(i) for i in eachcol(gmtx))
    mapreduce(hcat, eachrow(gmtx)) do smple
        _setcounts(smple, allalleles)
        #[j for i in _countset.(smple, allalleles) for j in i]
    end |> permutedims
end
precompile(countmatrix, (PopData,))

"""
    freqmatrix_zero(data::PopData)
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are replaced by zeros.
"""
function freqmatrix_zero(data::PopData)
    # divide each row (sample) by the ploidy of that sample
    out = countmatrix(data)
    replace!(out, -1 => 0)
    out ./ data.sampleinfo.ploidy
end
precompile(freqmatrix_zero, (PopData,))


"""
    freqmatrix_mean(data::PopData)
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are replaced by the global mean allele frequency.
"""
function freqmatrix_mean(data::PopData)
    counts = @inbounds countmatrix(data) ./ data.sampleinfo.ploidy
    map(eachcol(counts)) do alcol
        colmean = mean([x for x in alcol if x >= 0])
        replace!(x -> x < 0 ? colmean : x, alcol)
    end
    return counts
end
precompile(freqmatrix_mean, (PopData,))


"""
    freqmatrix_missing(data::PopData)
Create a matrix of allele frequencies per genotype where rows are samples
and columns are the frequency of an allele for that locus in that sample.
Missing values are kept as `missing`.
"""
function freqmatrix_missing(data::PopData)
    out = allowmissing(countmatrix(data))
    replace!(out, -1 => missing)
    out ./ data.sampleinfo.ploidy
end
precompile(freqmatrix_missing, (PopData,))


"""
    freqmatrix_scale(freqs::Matrix{Float32}, scale::Bool = true, center::Bool = true)
Returns a Z-score scaled matrix of allele frequencies where rows are samples 
and columns are the frequency of an allele for that locus in that sample.
- `scale`: a 'Bool' of whether to z-score scale allele frequencies (default: `false`)
- `center`: a `Bool` of whether to center the allele frequencies (default: 'false')
"""
function freqmatrix_scale(freqs::Matrix{Float64}, scale::Bool = true, center::Bool = true)
    mtx = standardize(ZScoreTransform, freqs, dims = 1, scale = scale, center = center)
    # replace almost-zero values caused by missing values with 0.0
    replace!(x ->  0 < x < (10^-9) ? 0.0 : x, mtx)
    return mtx
end
precompile(freqmatrix_scale, (PopData,))

"""
    featurematrix(data::PopData, matrixtype::String = "genotype")

### Positional Arguments
    - `data`: a PopData object
    - `matrixtype`: a `String` of `genotype`, or `allele` (default: `genotype`)

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
 julia> featurematrix(y, "allele")
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
function featurematrix(data::PopData, matrixtype::String = "genotype")::Matrix{Int8}
    if lowercase(matrixtype) == "genotype"
        featurematrix_genotype(data)
    elseif lowercase(matrixtype) == "allele"
        featurematrix_allele(data)
    else
        throw(ArgumentError("matrix type $matrixtype not recognized. Please use one of \"genotype\" or \"allele\"."))
    end
end

"""
    featurematrix_genotype(data::PopData, matrixtype::String = "genotype")
Return a matrix of dummy-encoded genotypes (0,1,2...), where rows correspond with samples and columns correspond to loci.
Missing genotypes are encoded as `-1`. For biallelic loci, `0` encodes homozygous for allele 1, `1` encodes for a heterozygote,
and `2` encodes for homozygous allele 2.
"""
function featurematrix_genotype(data::PopData)::Matrix{Int8}
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
precompile(freaturematrix_genotype, (PopData,))

"""
    featurematrix_allele(data::PopData)
Return a matrix of dummy-encoded alleles (0,1,2...), where rows correspond with samples and columns correspond to loci.
Missing genotypes are encoded as `-1`.
"""
function featurematrix_allele(data::PopData)::Matrix{Int8}
    allecounts = countmatrix(data)
    replace!(allecounts) do i
        i > 0 ? Int8(1) : i
    end
    return allecounts
end
precompile(freaturematrix_allele, (PopData,))