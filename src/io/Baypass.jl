"""
    _setcounts_biallelic(q, r)
Returns a vector of counts of alleles from Genotypes `q` in allele vector `r`.
This is distinct from `setcounts` in that `missing` values are preserved as `0` instead of `-1`
"""
function _setcounts_biallelic(q, r)
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
            @inbounds cnt[idx] = geno === missing ? 0 : count(==(h), geno)
        end
    end
    return cnt
end


"""
    countmatrix_biallelic(data::PopData)
Create a matrix of allele count per genotype where rows are samples
and columns are the occurence count of an allele for that locus in that sample.
`missing` values are preserved as `0``.
"""
function countmatrix_biallelic(data::PopData)
    gmtx = locimatrix(data)
    allalleles = Tuple(uniquealleles(i) for i in eachcol(gmtx))
    for i in allalleles
        # if there is only 1 allele, add -9 as a second (fake) allele
        # this will ensure two columns/alleles per locus in the output matrix
        length(i) == 1 ? push!(i, -9) : i
    end
    mapreduce(hcat, eachrow(gmtx)) do smple
        _setcounts_biallelic(smple, allalleles)
    end |> permutedims
end

"""
    baypass(data::PopData; filename::Union{String, Nothing} = nothing)
Convert a `PopData` object into a Baypass-format matrix. The required input format for the software
requires biallelic data. By default, it returns just the Baypass-format matrix; use the keyword argument `filename` to specify a file to write the matrix to.
This function **does not perform a Baypass analysis**, but instead creates the input matrix necessary for it.

The matrix specification is:
- rows = loci
    - each row is a different locus
- columns = allele counts per population
    - each pair of columns correspond to the alleles' counts (2 alleles, 2 columns) for a population
    - as a result, there should be 2 × n_populations columns
    - e.g. row 1, columns 1:2 are the allele counts for locus 1 in population 1

Baypass information: http://www1.montpellier.inra.fr/CBGP/software/baypass/
"""
function baypass(data::PopData; filename::Union{String, Nothing} = nothing)
    isbiallelic(data) || throw(ArgumentError("Like the Baypass software, the file writer requires biallelic data. Check that your data is biallelic at each locus or use dropmultiallelic!()"))
    pops = populations(data)
    bp_mtx = mapreduce(hcat, pops) do p
        # create new popdata by subsetting by column
        mtx = countmatrix_biallelic(data[data.genodata.population .== p])
        sums = sum(mtx, dims = 1)
        # get a list of the indices, split by odds and evens (odds = allele 1, evens = allele 2)
        idx = eachindex(sums)
        allele1 = sums[idx[isodd.(idx)]]
        allele2 = sums[idx[iseven.(idx)]]
        hcat(allele1, allele2)
    end
    if !isnothing(filename)
        if isfile(filename)
            @info "Overwriting $filename with Baypass data matrix"
            println()
        else
            @info "Writing Baypass data matrix to $filename"
            println()
        end
        open(filename, "w") do f
            for i in eachrow(bp_mtx)
                join(f, i, " ")
                println(f)
            end
        end
    end
    return bp_mtx
end