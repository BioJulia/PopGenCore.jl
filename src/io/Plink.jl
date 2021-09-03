
function plink(infile::String; silent::Bool = false)
    basefile = join(split(infile, ".")[1:end-1])
    famfile= CSV.read(
        basefile * ".fam", DataFrame, 
        header = [:population, :name, :sire, :dam, :sex, :phenotype],
        missingstrings = ["-9", "0"],
        types = Dict(:sex => Int8)
    )
    nsamples = length(famfile.name)
    npopulations = length(famfile.population)
    sirecheck = !all(ismissing.(famfile.sire))
    damcheck = !all(ismissing.(famfile.dam))
    sexcheck = !all(ismissing.(famfile.sex))
    phenotypecheck = !all(ismissing.(famfile.phenotype))
    # bimfile #
    # col1 chrom number
    # col2 snp name
    # col3 snp position in morgans/com
    # col4 snp physical position (resets for each chrom)
    # col5 and 6 the allele values usually 5=minor 6=major
    # 0 is missing
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

    CSV.read(basefile * ".bed", DataFrames)
end

