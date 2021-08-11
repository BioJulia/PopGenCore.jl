using .VariantCallFormat

export bcf, vcf

"""
    openvcf(::String)
Open VCF file (`.vcf(.gz)`, or `.bcf(.gz)`) and return an `IO` stream in reading mode `"r"`.
"""
function openvcf(infile::String)
    if endswith(infile, ".gz")
        throw(ArgumentError("Please load in GZip.jl with \`using GZip\`"))
    else
        return open(infile, "r")
    end
end

function bcf(infile::String; rename_loci::Bool = false, silent::Bool = false, allow_monomorphic::Bool = false)
    bases = (A = Int8(1), T = Int8(2), C = Int8(3), G = Int8(4), miss = Int8(0))
    stream = BCF.Reader(openvcf(infile))
    nmarkers = countlines(openvcf(infile)) - length(header(stream)) - 1
    sample_ID = header(stream).sampleID
    nsamples = length(sample_ID)
    geno_df = DataFrame(:name => sample_ID, :population =>  "missing")
    
    if silent == false
        @info "\n $(abspath(infile))\n data: samples = $nsamples, populations = 0, loci = $nmarkers\n ---> population info must be added"
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
    close(stream)
    if rename_loci
        rnm = vcat([:name, :population], [Symbol.("snp_" * i) for i in string.(1:nmarkers)])
        rename!(geno_df, rnm)
    end
    stacked_geno_df = DataFrames.stack(geno_df, DataFrames.Not(1:2))
    rename!(stacked_geno_df, [:name, :population, :locus, :genotype])
    # set columns as PooledArrays
    select!(
        stacked_geno_df, 
        :name => (i -> PooledArray(i, compress = true)) => :name, 
        :population => (i -> PooledArray(i, compress = true)) => :population, 
        :locus => (i -> PooledArray(i |> Vector{String}, compress = true)) => :locus, 
        :genotype
    )
    # replace missing genotypes as missing
    stacked_geno_df.genotype = map(i -> all(0 .== i) ? missing : i, stacked_geno_df.genotype)
    sort!(stacked_geno_df, [:name, :locus])
    meta_df = generate_meta(stacked_geno_df)
    pd_out = PopData(meta_df, stacked_geno_df)
    pd_out = !allow_monomorphic ? drop_monomorphic(pd_out, silent = silent) : pd_out
    return pd_out
end

### VCF parsing ###

function vcf(infile::String; rename_loci::Bool = false, silent::Bool = false, allow_monomorphic::Bool = false)
    bases = (A = Int8(1), T = Int8(2), C = Int8(3), G = Int8(4), miss = Int8(0))
    stream = VCF.Reader(openvcf(infile))
    nmarkers = countlines(openvcf(infile)) - length(header(stream)) - 1
    sample_ID = header(stream).sampleID
    nsamples = length(sample_ID)
    geno_df = DataFrame(:name => sample_ID, :population =>  "missing")
    
    if silent == false
        @info "\n $(abspath(infile))\n data: samples = $nsamples, populations = 0, loci = $nmarkers\n ---> population info must be added"
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
    close(stream)
    if rename_loci
        rnm = vcat([:name, :population], [Symbol.("snp_" * i) for i in string.(1:nmarkers)])
        rename!(geno_df, rnm)
    end
    stacked_geno_df = DataFrames.stack(geno_df, DataFrames.Not(1:2))
    rename!(stacked_geno_df, [:name, :population, :locus, :genotype])
    # set columns as PooledArrays
    select!(
        stacked_geno_df, 
        :name => (i -> PooledArray(i, compress = true)) => :name, 
        :population => (i -> PooledArray(i, compress = true)) => :population, 
        :locus => (i -> PooledArray(i |> Vector{String}, compress = true)) => :locus, 
        :genotype
    )
    # replace missing genotypes as missing
    stacked_geno_df.genotype = map(i -> all(0 .== i) ? missing : i, stacked_geno_df.genotype)
    sort!(stacked_geno_df, [:name, :locus])
    meta_df = generate_meta(stacked_geno_df)
    pd_out = PopData(meta_df, stacked_geno_df)
    pd_out = !allow_monomorphic ? drop_monomorphic(pd_out, silent = silent) : pd_out
    return pd_out
end

#vcf(normpath(joinpath(@__DIR__, "precompile", "precompile")) * ".vcf", silent = true) ;