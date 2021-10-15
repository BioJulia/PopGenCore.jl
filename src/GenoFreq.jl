export geno_freq, geno_freq_expected, geno_count_observed, geno_count_expected

"""
    geno_count_observed(locus::GenoArray)
Return a `Dict` of genotype counts of a single locus in a
`PopData` object.
"""
@inline function geno_count_observed(locus::T) where T<:GenoArray
    # conditional testing if all genos are missing
    all(ismissing.(locus)) && return missing
    countmap(skipmissing(locus))
end

"""
    geno_count_expected(locus::GenoArray)
Return a `Dict` of the expected genotype counts of a single locus in a
`PopData` object. Expected counts are calculated as the product of observed
allele frequencies multiplied by the number of non-missing genotypes.
"""
function geno_count_expected(locus::T) where T<:GenoArray
    #count number of non-missing genotypes in the locus
    n = nonmissing(locus)

    # Get expected number of genotypes in a locus
    ## get the observed allele frequencies
    allele_dict = allele_freq(locus)

    ## split the appropriate pairs into their own vectors
    alle, freq = collect(keys(allele_dict)), collect(values(allele_dict))
    ## calculate expected genotype frequencies by multiplying all-by-all x n
    expected_genotype_freq = vec(freq * freq' .* n)

    # reform genotype frequencies with same all-by-all approach
    genos = reverse.(Base.Iterators.product(alle, alle) |> collect |> vec)
    expected = Dict{nonmissingtype(eltype(locus)), Float64}()
    for (geno, freq) in zip(genos, expected_genotype_freq)
        expected[geno] = get!(expected, geno, 0.0) + freq
    end

    return expected
end


"""
    geno_freq(locus::GenoArray)
Return a `Dict` of genotype frequencies of a single locus in a
`PopData` object.
"""
@inline function geno_freq(locus::T) where T<:GenoArray
    # conditional testing if all genos are missing
    all(ismissing.(locus)) && return missing
    proportionmap(locus |> skipmissing |> collect)
end


"""
    geno_freq(data::PopData, locus::String; population::Bool = false)
Return a `Dict` of genotype frequencies of a single locus in a `PopData`
object. Use `population = true` to return a table of genotype frequencies
of that locus per population.
### Example
```
cats = @nancycats;
geno_freq(cats, "fca8")
geno_freq(cats, "fca8", population = true)
```
"""
function geno_freq(data::PopData, locus::String; population::Bool=false)
    if !population
        tmp = data.genodata[data.genodata.locus .== locus, :genotype]
        geno_freq(tmp)
    else
        tmp = data.genodata[data.genodata.locus .== locus, :]
        DataFrames.combine(groupby(tmp, :population), :genotype => geno_freq => :freq)
    end
end

"""
    geno_freq_expected(locus::GenoArray)
Return a `Dict` of the expected genotype frequencies of a single locus in a
`PopData` object. Expected frequencies are calculated as the product of
observed allele frequencies.
"""
function geno_freq_expected(locus::T) where T<:GenoArray
    # Get expected number of genotypes in a locus
    ## get the observed allele frequencies
    allele_dict = allele_freq(locus)

    ## split the appropriate pairs into their own vectors
    alle, freq = collect(keys(allele_dict)), collect(values(allele_dict))

    ## calculate expected genotype frequencies by multiplying all-by-all
    expected_genotype_freq = vec(freq * freq')

    # reform genotype frequencies with same all-by-all approach
    genos = reverse.(Base.Iterators.product(alle, alle) |> collect |> vec)

    # reform genotype frequencies into a Dict
    expected = Dict{nonmissingtype(eltype(locus)), Float64}()
    for (geno, freq) in zip(genos, expected_genotype_freq)
        expected[geno] = get!(expected, geno, 0.0) + freq
    end

    return expected
end

"""
    geno_freq_expected(data::PopData, locus::String; population::Bool = false)
Return a `Dict` of expected genotype frequencies of a single locus in a
`PopData` object. Use `population = true` to return a table of expected genotype
frequencies of that locus per population.
### Example
```
cats = @nancycats;
geno_freq_expected(cats, "fca8")
geno_freq_expected(cats, "fca8", population = true)
```
"""
function geno_freq_expected(data::PopData, locus::String; population::Bool=false)
    if !population
        tmp = data.genodata[data.genodata.locus .== locus, :genotype]
        geno_freq_expected(tmp)
    else
        tmp = data.genodata[data.genodata.locus .== locus, :]
        DataFrames.combine(groupby(tmp, :population), :genotype => geno_freq_expected => :freq)
    end
end
