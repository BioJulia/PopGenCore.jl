export bcf, vcf

"""
    openvcf(::String)
Open Variant Call Format file `.vcf[.gz]`, or `.bcf[.gz]` and return an `IO` stream in reading mode `"r"`.
"""
function openvcf(infile::String)
    if endswith(infile, ".gz")
        return GzipDecompressorStream(open(infile, "r"))
    else
        return open(infile, "r")
    end
end

"""
    bcf(infile::String; ; rename_loci::Bool, silent::Bool, allow_monomorphic::Bool)
Load a BCF file into memory as a PopData object. Population information needs to be provided separately. 
- `infile` : path to BCF file (can be gzipped)

**Keyword Arguments**
- `rename_loci` : true/false of whether to simplify loci names to "snp_#" (default: `false`)
- `allow_monomorphic` : true/false of whether to keep monomorphic loci (default: `false`)
- `silent`: true/false of whether to print extra file information (default: `false`).
Alleles are recoded according to the following schema:


| **Base**   |  A   |  T   |  C   |  G   |
| :--------  | :--: | :--: | :--: | :--: |
| **Allele** |  1   |  2   |  3   |  4   |


### Mixed-ploidy data
If importing mixed-ploidy data (such as poolseq), you will need to perform an additional
step to convert the genotype column into the correct `GenoArray` type:
```julia
julia> mydata = bcf("path/to/file.bcf", silent = true, rename_loci = true) ;

julia> mydata.genodata.genotype =  mydata.genodata.genotype |> Array{Union{Missing, NTuple}}
```
"""
function bcf(infile::String; rename_loci::Bool = false, silent::Bool = false, allow_monomorphic::Bool = false)
    bases = (A = Int8(1), T = Int8(2), C = Int8(3), G = Int8(4), miss = Int8(0))
    counts = countlines(openvcf(infile))
    f = openvcf(infile)
    stream = BCF.Reader(f)
    nmarkers = counts - length(header(stream)) - 1
    sample_ID = header(stream).sampleID
    nsamples = length(sample_ID)
    geno_df = DataFrame(:name => sample_ID, :population =>  "missing")
    if silent == false
        @info "\n $(abspath(infile))\n data: samples = $nsamples, populations = 0, loci = $nmarkers\n ⟶  population info must be added"
        println()
    end

    for record in stream
        ref_alt = Dict(-1 => "miss", 0 => BCF.ref(record))
        [ref_alt[i] = j for (i,j) in enumerate(BCF.alt(record))]
        raw_geno = BCF.genotype(record, 1:nsamples, "GT")
        conv_geno = map(raw_geno) do rg
            tmp = replace.(rg, "." => "-1")
            ig = collect(parse.(Int8, split(tmp, r"\/|\|")))
            [bases[Symbol(ref_alt[i])] for i in ig] |> sort |> Tuple
        end
        insertcols!(geno_df, Symbol(BCF.chrom(record) * "_" * string(BCF.pos(record))) => conv_geno)
    end
    close(stream) ; close(f)

    if rename_loci
        rnm = vcat([:name, :population], [Symbol.("snp_" * i) for i in string.(1:nmarkers)])
        rename!(geno_df, rnm)
    end

    geno_df = DataFrames.stack(geno_df, DataFrames.Not(1:2))
    rename!(geno_df, [:name, :population, :locus, :genotype])
    
    # set columns as PooledArrays
    select!(
        geno_df, 
        :name => (i -> PooledArray(i, compress = true)) => :name, 
        :population => (i -> PooledArray(i, compress = true)) => :population, 
        :locus => (i -> PooledArray(i |> Vector{String}, compress = true)) => :locus, 
        :genotype => (j -> map(i -> all(0 .== i) ? missing : i, j)) => :genotype
    )

    # replace missing genotypes as missing
    sort!(geno_df, [:locus, :population, :name], lt = natural)
    meta_df = generate_meta(geno_df)
    pd_out = PopData(meta_df, geno_df)
    !allow_monomorphic && drop_monomorphic!(pd_out, silent = silent)
    return pd_out
end

### VCF parsing ###


"""
    vcf(infile::String; ; rename_loci::Bool, silent::Bool, allow_monomorphic::Bool)
Load a VCF file into memory as a PopData object. Population information needs to be provided separately. 
- `infile` : path to VCF file (can be gzipped)

**Keyword Arguments**
- `rename_loci` : true/false of whether to simplify loci names to "snp_#" (default: `false`)
- `allow_monomorphic` : true/false of whether to keep monomorphic loci (default: `false`)
- `silent`: true/false of whether to print extra file information (default: `false`).
Alleles are recoded according to the following schema:


| **Base**   |  A   |  T   |  C   |  G   |
| :--------  | :--: | :--: | :--: | :--: |
| **Allele** |  1   |  2   |  3   |  4   |


### Mixed-ploidy data
If importing mixed-ploidy data (such as poolseq), you will need to perform an additional
step to convert the genotype column into the correct `GenoArray` type:
```julia
julia> mydata = vcf("path/to/file.vcf", silent = true, rename_loci = true) ;

julia> mydata.genodata.genotype =  mydata.genodata.genotype |> Array{Union{Missing, NTuple}}

```
"""
function vcf(infile::String; rename_loci::Bool = false, silent::Bool = false, allow_monomorphic::Bool = false)
    bases = (A = Int8(1), T = Int8(2), C = Int8(3), G = Int8(4), miss = Int8(0))
    counts = countlines(openvcf(infile))
    f = openvcf(infile)
    stream = VCF.Reader(f)
    nmarkers = counts - length(header(stream)) - 1
    sample_ID = header(stream).sampleID
    nsamples = length(sample_ID)
    geno_df = DataFrame(:name => sample_ID, :population =>  "missing")
    if silent == false
        @info "\n $(abspath(infile))\n data: samples = $nsamples, populations = 0, loci = $nmarkers\n ⟶  population info must be added"
        println()
    end

    for record in stream
        ref_alt = Dict(-1 => "miss", 0 => VCF.ref(record))
        [ref_alt[i] = j for (i,j) in enumerate(VCF.alt(record))]
        raw_geno = VCF.genotype(record, 1:nsamples, "GT")
        conv_geno = map(raw_geno) do rg
            tmp = replace.(rg, "." => "-1")
            ig = collect(parse.(Int8, split(tmp, r"\/|\|")))
            [bases[Symbol(ref_alt[i])] for i in ig] |> sort |> Tuple
        end
        insertcols!(geno_df, Symbol(VCF.chrom(record) * "_" * string(VCF.pos(record))) => conv_geno)
    end
    close(stream) ; close(f)

    if rename_loci
        rnm = vcat([:name, :population], [Symbol.("snp_" * i) for i in string.(1:nmarkers)])
        rename!(geno_df, rnm)
    end

    geno_df = DataFrames.stack(geno_df, DataFrames.Not(1:2))
    rename!(geno_df, [:name, :population, :locus, :genotype])
    
    # set columns as PooledArrays
    select!(
        geno_df, 
        :name => (i -> PooledArray(i, compress = true)) => :name, 
        :population => (i -> PooledArray(i, compress = true)) => :population, 
        :locus => (i -> PooledArray(i |> Vector{String}, compress = true)) => :locus, 
        :genotype => (j -> map(i -> all(0 .== i) ? missing : i, j)) => :genotype
    )

    # replace missing genotypes as missing
    sort!(geno_df, [:locus, :population, :name], lt = natural)
    meta_df = generate_meta(geno_df)
    pd_out = PopData(meta_df, geno_df)
    !allow_monomorphic && drop_monomorphic!(pd_out, silent = silent)
    return pd_out
end