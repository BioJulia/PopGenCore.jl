
function plink(infile::String; silent::Bool = false)
    basefile = join(split(infile, ".")[1:end-1])
    famfile = CSV.read(
        basefile * ".fam", DataFrame, 
        header = [:population, :name, :sire, :dam, :sex, :phenotype],
        missingstrings = ["-9", "0"],
        types = Dict(:sex => Int8)
    )
    nsamples = length(famfile.name)
    npopulations = length(unique(famfile.population))
    sirecheck = !all(ismissing.(famfile.sire))
    damcheck = !all(ismissing.(famfile.dam))
    sexcheck = !all(ismissing.(famfile.sex))
    phenotypecheck = !all(ismissing.(famfile.phenotype))
    bimfile = CSV.read(
        basefile * ".bim",
        DataFrame,
        header = [:chrom, :snp, :cM, :bp, :allele1, :allele2],
        missingstrings = ["0"],
        drop = [3,4],
        types = Dict(:chrom => String, :allele1 => Char, :allele2 => Char),
        pool=0.3
    )
    locinames = bimfile.chrom .* "_" .* bimfile.snp
    nloci =  length(locinames)

    if !silent
        @info "\n $(abspath(infile))\n data: loci = $(nloci), samples = $(nsamples), populations = $(npopulations)\n other: sire = $(sirecheck), dam = $(damcheck), sex = $(sexcheck), phenotype = $(phenotypecheck)"
        println()
    end

    bedfile = basefile * ".bed"
    data = open(bedfile, "r") do io
        read(io, UInt16) == 0x1b6c || throw(ArgumentError("wrong magic number in file $bedfile"))
        read(io, UInt8) == 0x01 || throw(ArgumentError(".bed file, $bedfile, is not in correct orientation"))
        if endswith(bedfile, ".bed")
            return Mmap.mmap(io)
        else
            return read(io)
        end
    end
    nrows = (nsamples + 3) >> 2   # the number of rows in the Matrix{UInt8}
    n, r = divrem(length(data), nrows)
    iszero(r) || throw(ArgumentError("filesize of $bedfile is not a multiple of $nrows"))
    reshape(data, (nrows, n))
    return famfile
end