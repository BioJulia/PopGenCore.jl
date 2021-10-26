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


function _plinkped(infile::String, keepfields::Union{Symbol,Vector{Symbol}} = :all, silent::Bool = false)
    basefile = splitext(infile)[1]
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
        locinames = bimfile.snp
    end
    pedfile = CSV.read(
        infile, 
        DataFrame,  
        header = false,
        missingstring = ["-9", "0"],
        typemap = Dict(Int64 => Int8)
        )
    n = size(pedfile)[2] - 6
    locinames = bimfound ? locinames : ["snp_" * "$i" for i in 1:n]
    rename!(pedfile, vcat([:population, :name, :sire, :dam, :sex, :phenotype], Symbol.(locinames)))
    return pedfile
    nsamples = length(pedfile.name)
    npopulations = length(unique(pedfile.population))
    sirecheck = !all(ismissing.(pedfile.sire)) ?  "(✔)" : "(𐄂)"
    damcheck = !all(ismissing.(pedfile.dam)) ?  "(✔)" : "(𐄂)"
    sexcheck = !all(ismissing.(pedfile.sex)) ?  "(✔)" : "(𐄂)"
    phenotypecheck = !all(ismissing.(pedfile.phenotype)) ?  "(✔)" : "(𐄂)"
    if !silent
        if !bimfound
            @info "\n $(truncatepath(abspath(infile)))\n data: loci = $(n), samples = $(nsamples), populations = $(npopulations)\n .fam: sire $(sirecheck), dam $(damcheck), sex $(sexcheck), phenotype $(phenotypecheck)\n $(basefile * ".bim") not found, generating new marker names"
        else
            @info "\n $(truncatepath(abspath(infile)))\n data: loci = $(n), samples = $(nsamples), populations = $(npopulations)\n .fam: sire $(sirecheck), dam $(damcheck), sex $(sexcheck), phenotype $(phenotypecheck)"
        end
        println()
    end
end


function _plinkbed(infile::String, famfields::Union{Symbol,Vector{Symbol}} = :all, bimfields::Union{Symbol,Vector{Symbol}} = :all, silent::Bool = false)
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
            header = [:chromosome, :snp, :cm, :bp, :allele1, :allele2],
            missingstring = ["0"],
            #drop = [3,4],
            types = Dict(:chromosome => String, :allele1 => Char, :allele2 => Char),
            pool=0.3
        )
        locinames = bimfile.snp
        chromcheck = !all(ismissing.(bimfile.chromosome)) ?  "(✔)" : "(𐄂)"
        cmcheck = !all(ismissing.(bimfile.cm)) ?  "(✔)" : "(𐄂)"
        bpcheck = !all(ismissing.(bimfile.bp)) ?  "(✔)" : "(𐄂)"
    end
    bedfile = basefile * ".bed"
    data = open(bedfile, "r") do io
        Base.read(io, UInt16) == 0x1b6c || throw(ArgumentError("Incorrect \"magic\" number in $bedfile. It should be 0x1b6c"))
        Base.read(io, UInt8) == 0x01 || throw(ArgumentError("$bedfile is not in the correct orientation"))
        return Base.read(io)
    end
    nrows = (nsamples + 3) >> 2   # the number of rows in the Matrix{UInt8}
    n, r = divrem(length(data), nrows)
    locinames = bimfound ? locinames : ["snp_" * "$i" for i in 1:n]
    iszero(r) || throw(ArgumentError("The filesize of $bedfile is not a multiple of $nrows (and it should be)"))
    if !silent
        if !bimfound
            @info "\n $(truncatepath(abspath(infile)))\n data: loci = $(n), samples = $(nsamples), populations = $(npopulations)\n .fam: sire $(sirecheck), dam $(damcheck), sex $(sexcheck), phenotype $(phenotypecheck)\n $(basefile * ".bim") not found, generating new marker names"
        else
            @info "\n $(truncatepath(abspath(infile)))\n data: loci = $(n), samples = $(nsamples), populations = $(npopulations)\n .fam: sire $(sirecheck), dam $(damcheck), sex $(sexcheck), phenotype $(phenotypecheck)\n .bim: chromosome $(chromcheck), position $(cmcheck), bp $(bpcheck)"
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
    #TODO generate_meta is deprecated
    metadf = generate_meta(genodf)
    if famfields != :none
        if famfields == :all
            sirecheck == "(✔)" && insertcols!(metadf, :sire => famfile.sire)
            damcheck == "(✔)" && insertcols!(metadf, :dam => famfile.dam)
            sexcheck == "(✔)" && insertcols!(metadf, :sex => famfile.sex)
            phenotypecheck == "(✔)" && insertcols!(metadf, :phenotype => famfile.phenotype)
        else (typeof(famfields) == Symbol) | (typeof(famfields) == Vector{Symbol})
            _fields = typeof(famfields) == Symbol ? [famfields] : famfields
            checkfields = symdiff(_fields, intersect(_fields, [:sire, :dam, :sex, :phenotype]))
            length(checkfields) > 0 && throw(ArgumentError("Unknown .fam file field(s) requested: $checkfields"))
            [insertcols!(metadf, i => famfile[:, i]) for i in _fields]
        end
    end
    if bimfound && bimfields != :none
        if bimfields == :all
            chromcheck == "(✔)" && insertcols!(genodf, 3, :chromosome => PooledArray(repeat(coalesce.(bimfile.chromosome,0), inner = nsamples), compress = true))
            cmcheck == "(✔)" && insertcols!(genodf, :cm => PooledArray(repeat(Float64.(coalesce.(bimfile.cm, 0)), inner = nsamples), compress = true))
            bpcheck == "(✔)" &&  insertcols!(genodf, :bp => PooledArray(repeat(Int64.(coalesce.(bimfile.bp, 0)), inner = nsamples), compress = true))
        else (typeof(bimfields) == Symbol) | (typeof(bimfields) == Vector{Symbol})
            _fields = typeof(bimfields) == Symbol ? [bimfields] : bimfields
            checkfields = symdiff(_fields, intersect(_fields, [:chromosome, :cm, :bp]))
            length(checkfields) > 0 && throw(ArgumentError("Unknown .bim file field(s) requested: $checkfields"))
            :chromosome ∈ _fields && insertcols!(genodf, 3, :chromosome => PooledArray(repeat(coalesce.(bimfile.chromosome,0), inner = nsamples), compress = true))
            :cm ∈ _fields && insertcols!(genodf, :cm => PooledArray(repeat(Float64.(coalesce.(bimfile.cm, 0)), inner = nsamples), compress = true))
            :bp ∈ _fields &&  insertcols!(genodf, :bp => PooledArray(repeat(Int64.(coalesce.(bimfile.bp, 0)), inner = nsamples), compress = true))
            #[insertcols!(genodf, i => PooledArray(repeat(coalesce.(bimfile[:, i], 0), inner = nsamples), compress = true)) for i in _fields]
        end
    end
    return PopData(metadf, genodf)
end

"""
    plink(infile::String; keepfields::Symbol|Vector{Symbol}, silent::Bool)
Read a PLINK `.ped` or binary `.bed` file into memory as a `PopData` object.
Requires an accompanying `.fam` file in the same directory, but an accompanying `.bim` file is optional.
- `infile::String` : path to `.ped` or `.bed` file

### Keyword Arguments
- `famfields::Symbol|Vector{Symbol}` : which additional fields to import from the `.fam` file
    - `:all` [default]
    - `:none`
    - any one or combination of `[:sire, :dam, :sex, :phenotype]`
- `bimfields::Symbol|Vector{Symbol}` : which additional fields to import from the optional `.bim` file
    - `:all` [default]
    - `:none`
    - any one or combination of `[:chromosome, :cm, :bp]`
- `silent::Bool`   : whether to print file information during import (default: `false`)

## Example
```julia
para = plink("datadir/parakeet.ped", famfields = :sex)
parr = plink("datadir/parrot.bed", famfields = [:sire, :dam], bimfields = :chromosome)
```
"""
function plink(infile::String; famfields::Union{Symbol,Vector{Symbol}} = :all, bimfields::Union{Symbol,Vector{Symbol}} = :all, silent::Bool = false)
    isfile(infile) || throw(ArgumentError("$infile not found in working directory."))
    if endswith(infile, ".ped")
        _plinkped(infile, famfields, bimfields, silent)
    elseif endswith(infile, ".bed")
        _plinkbed(infile, famfields, bimfields, silent)
    else
        throw(ArgumentError("Filename $infile is not recognized as having a .ped or .bed extension"))
    end
end

### writing to PLINK format ###

@inline _genoconversion(genotype::T) where T<:Genotype = join(genotype, " ")
@inline _genoconversion(genotype::Missing) = "0 0"


"""
    plink(data::PopData; filename::String)

Write a biallelic `PopData` object to PLINK `.ped` format with an accompanying
`.fam` file. Genotypes are coded by the PLINK standard:
- Integers are the alleles
- `0` encodes missing
- After column 6, every two numbers indicate a diploid genotype.
## Example
```julia
sharks = dropmultiallelic(@gulfsharks) ;
plink(sharks, filename = "biallelic_sharks.ped")
```
"""
function plink(data::PopData; filename::String)
    data.metadata.biallelic != true && throw(ArgumentError("To write to PLINK format, data must be biallelic.\nThis can be done with dropmultiallelic(::PopData)\n"))
    basename = endswith(filename, r".bed|.ped") ? splitext(filename)[1] : filename
    # the .fam file
    tmp = data.metadata[:, r"population|name|sire|dam|sex|phenotype"]
    current_cols = propertynames(tmp)
    fill_cols = symdiff(current_cols,[:population, :name, :sire, :dam, :sex, :phenotype])
    for i in fill_cols
        tmp[!, i] .= Int8(0)
    end
    select!(tmp, :population, :name, :sire, :dam, :sex, :phenotype)
    tmp.sex = coalesce.(tmp.sex, 0)
    CSV.write(basename * ".fam", tmp, delim = " ", writeheader = false, quotestrings = false, missingstring = "0")
    # the .ped file
    open(filename, "w") do outped 
        eachind = groupby(data.genodata, :name)
        eachindmeta = groupby(tmp, :name)
        for i in data.genodata.name.pool
            metarow = eachindmeta[(;name = i)]
            metainfo = [metarow.population[1], metarow.name[1], metarow.sire[1], metarow.dam[1], metarow.sex[1], metarow.phenotype[1]]
            genos = eachind[(;name = i)].genotype
            printrow = vcat(metainfo, _genoconversion.(genos))
            join(outped, printrow, " ")
            println(outped, "")
        end
    end
    # the bim file
    tmp = unique(data.genodata, :locus)[:, r"chrom|locus|cm|bp"]
    current_cols = propertynames(tmp)
    fill_cols = symdiff(current_cols,[:chrom, :locus, :cm, :bp])
    for i in fill_cols
        tmp[!, i] .= Int8(0)
    end
    select!(tmp, :chrom, :locus, :cm, :bp)
    CSV.write(basename * ".bim", tmp, delim = " ", writeheader = false, quotestrings = false)
    return
end


# BED FORMAT
"""
    plink(data::PopData; filename::String)

Write a biallelic `PopData` object to PLINK `.bed` format with an accompanying
`.fam` file. Genotypes are coded by the PLINK standard:
- `00` Homozygous for first allele (`0x00`)
- `01` Missing genotype (`0x01`)
- `10` Heterozygous  (`0x02`)
- `11` Homozygous for second allele (`0x03`)

## Example
```julia
sharks = dropmultiallelic(@gulfsharks) ;
plink(sharks, filename = "biallelic_sharks.ped")
```
"""