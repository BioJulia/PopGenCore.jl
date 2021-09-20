# this function copies the getindex tool OpenMendel/SnpArrays.jl uses 
# to pull out the byte values from the compressed hex genotypes
# Repo: https://github.com/OpenMendel/SnpArrays.jl
@inline function _plinkindex(s::Matrix{UInt8}, i::Integer, j::Integer)
    #@boundscheck checkbounds(s, i, j)
    ip3 = i + 3
    (s[ip3 >> 2, j] >> ((ip3 & 0x03) << 1)) & 0x03
end


@inline function _plinkindex(s::Matrix{UInt8})
    _dim1, _dim2 = size(s)
    @inbounds reshape([_plinkindex(s, i, j) for i in 1:(_dim1*4) for j in 1:_dim2], :, _dim2)
end

@inline function _SNP(genotype::UInt8)
    #=
    00	Homozygous for first allele (0x00)
    01	Missing genotype (0x01)
    10	Heterozygous  (0x02)
    11	Homozygous for second allele in .bim file (0x03)
    =#
    genotype == 0x00 ? NTuple{2,Int8}((Int8(1),Int8(1))) :
        genotype == 0x01 ? missing :
        genotype == 0x02 ? NTuple{2,Int8}((Int8(1),Int8(2))) : NTuple{2,Int8}((Int8(2),Int8(2))) 
end

@inline function _SNP(genomatrix::AbstractArray{UInt8})
    Matrix{Union{Missing, NTuple{2,Int8}}}(_SNP.(genomatrix))
end


"""
    plink(infile::String; keepfields::Symbol|Vector{Symbol}, silent::Bool)
Read a PLINK `.ped` or binary `.bed` file into memory as a `PopData` object.
Requires an accompanying `.fam` file in the same directory, but an accompanying `.bim` file is optional.
- `infile::String` : path to `.ped` or `.bed` file

### Keyword Arguments
- `keepfields::Symbol|Vector{Symbol}` : which additional fields to import from the `.fam` file)
    - `:all` [default]
    - `:none`
    - any one or combination of `[:sire, :dam, :sex, :phenotype]`
- `silent::Bool`   : whether to print file information during import (default: `false`)

## Example
```julia
para = plink("datadir/parakeet.ped", keepfields = :sex)
parr = plink("datadir/parrot.bed", keepfields = [:sire, :dam])
```
"""
function plink(infile::String; keepfields::Union{Symbol,Vector{Symbol}} = :all, silent::Bool = false)
    isfile(infile) || throw(ArgumentError("$infile not found."))
    basefile = splitext(infile)[1]
    !isfile(basefile * ".fam") && throw(ArgumentError("$(basefile * ".fam") is required but wasn't found")) 
    famfile = CSV.read(
        basefile * ".fam", DataFrame, 
        header = [:population, :name, :sire, :dam, :sex, :phenotype],
        missingstring = ["-9", "0"],
        types = Dict(:sex => Int8)
    )
    nsamples = length(famfile.name)
    npopulations = length(unique(famfile.population))
    sirecheck = !all(ismissing.(famfile.sire)) ?  "(✔)" : "(𐄂)"
    damcheck = !all(ismissing.(famfile.dam)) ?  "(✔)" : "(𐄂)"
    sexcheck = !all(ismissing.(famfile.sex)) ?  "(✔)" : "(𐄂)"
    phenotypecheck = !all(ismissing.(famfile.phenotype)) ?  "(✔)" : "(𐄂)"
    bimfound = isfile(basefile * ".bim")
    if bimfound
        bimfile = CSV.read(
            basefile * ".bim",
            DataFrame,
            header = [:chrom, :snp, :cM, :bp, :allele1, :allele2],
            missingstring = ["0"],
            drop = [3,4],
            types = Dict(:chrom => String, :allele1 => Char, :allele2 => Char),
            pool=0.3
        )
        locinames = bimfile.chrom .* "_" .* bimfile.snp
    end
    if endswith(infile, ".bed")
        bedfile = basefile * ".bed"
        data = open(bedfile, "r") do io
            read(io, UInt16) == 0x1b6c || throw(ArgumentError("Incorrect \"magic\" number in $bedfile. It should be 0x1b6c"))
            read(io, UInt8) == 0x01 || throw(ArgumentError("$bedfile is not in the correct orientation"))
            return read(io)
        end
        nrows = (nsamples + 3) >> 2   # the number of rows in the Matrix{UInt8}
        n, r = divrem(length(data), nrows)
        iszero(r) || throw(ArgumentError("The filesize of $bedfile is not a multiple of $nrows (and it should be)"))
        if !silent
            if !bimfound
                @info "\n $(truncatepath(abspath(infile)))\n data: loci = $(n), samples = $(nsamples), populations = $(npopulations)\n .fam: sire $(sirecheck), dam $(damcheck), sex $(sexcheck), phenotype $(phenotypecheck)\n $(basefile * ".bim") not found, generating new marker names"
            else
                @info "\n $(truncatepath(abspath(infile)))\n data: loci = $(n), samples = $(nsamples), populations = $(npopulations)\n .fam: sire $(sirecheck), dam $(damcheck), sex $(sexcheck), phenotype $(phenotypecheck)"
            end
            println()
        end
        locinames = bimfound ? locinames : ["snp_" * "$i" for i in 1:n]
        genodf = DataFrame(
            :name => PooledArray(repeat(famfile.name, n) , compress = true),
            :population => PooledArray(repeat(famfile.population, n) , compress = true),
            :locus => PooledArray(repeat(locinames, inner = nsamples) , compress = true),
            :genotype => vec(_SNP(_plinkindex(reshape(data, (nrows, n))))) 
        )
    else
        data = CSV.read(basefile * ".ped", DataFrame, header = false, drop =1:6)
        # TODO this is incomplete
    end
    metadf = generate_meta(genodf)
    if keepfields != :none
        if keepfields == :all
            insertcols!(metadf, :sex => famfile.sex, :dam => famfile.dam, :sire => famfile.sire, :phenotype => famfile.phenotype)
        else (typeof(keepfields) == Symbol) | (typeof(keepfields) == Vector{Symbol})
            _fields = typeof(keepfields) == Symbol ? [keepfields] : keepfields
            checkfields = symdiff(_fields, intersect(_fields, [:sire, :dam, :sex, :phenotype]))
            length(checkfields) > 0 && throw(ArgumentError("Unknown .fam file field(s) requested: $checkfields"))
            [insertcols!(metadf, i => famfile[:, i]) for i in _fields]
        end
    end
    return PopData(metadf, genodf)
end

### writing to PLINK format ###

@inline function _genoconversion(genotype::T) where T<:Genotype
    #=
    00	Homozygous for first allele (0x00)
    01	Missing genotype (0x01)
    10	Heterozygous  (0x02)
    11	Homozygous for second allele in .bim file (0x03)
    =#
    genotype == NTuple{2,Int8}((Int8(1),Int8(1))) ? 0x00 :
        genotype == NTuple{2,Int8}((Int8(1),Int8(2))) ? 0x02 : 0x03
end

@inline function _genoconversion(genotype::Missing)
    return 0x01
end

"""
    plink(data::PopData; filename::String)

Write a biallelic `PopData` object to PLINK `.ped` format with an accompanying
`.fam` file. Genotypes are coded by the PLINK standard:
- `00` Homozygous for first allele (`0x00`)
- `01` Missing genotype (`0x01`)
- `10` Heterozygous  (`0x02`)
- `11` Homozygous for second allele (`0x03`)

## Example
```julia
sharks = drop_multiallelic(@gulfsharks) ;
plink(sharks, filename = "biallelic_sharks.ped")
```
"""
function plink(data::PopData; filename::String)
    data.info.biallelic != true && throw(ArgumentError("To write to PLINK format, data must be biallelic.\nThis can be done with drop_multiallelic(::PopData)\n"))
    basename = endswith(filename, r".bed|.ped") ? splitext(filename)[1] : filename
    # the .fam file
    tmp = data.metadata[:, r"population|name|sire|dam|sex|phenotype"]
    current_cols = propertynames(tmp)
    fill_cols = symdiff(current_cols,[:population, :name, :sire, :dam, :sex, :phenotype])
    for i in fill_cols
        tmp[!, i] .= Int8(0)
    end
    select!(tmp, :population, :name, :sire, :dam, :sex, :phenotype)
    CSV.write(basename * ".fam", tmp, delim = " ", writeheader = false, quotestrings = false)
    # the .ped file
    genomtx = loci_matrix(data)
    converted_geno = mapreduce(hcat, eachcol(genomtx)) do col
        _genoconversion.(col)
    end
    ped_df = hcat(tmp, DataFrame(converted_geno, :auto))
    CSV.write(basename * ".ped", ped_df, delim = " ", writeheader = false, quotestrings = false)
    return
end