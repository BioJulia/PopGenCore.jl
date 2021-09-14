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
                @info "\n $(abspath(infile))\n data: loci = $(n), samples = $(nsamples), populations = $(npopulations)\n .fam: sire $(sirecheck), dam $(damcheck), sex $(sexcheck), phenotype $(phenotypecheck)\n $(basefile * ".bim") not found, generating new marker names"
            else
                @info "\n $(abspath(infile))\n data: loci = $(n), samples = $(nsamples), populations = $(npopulations)\n .fam: sire $(sirecheck), dam $(damcheck), sex $(sexcheck), phenotype $(phenotypecheck)"
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
        data = CSV.read(infile * ".ped", header = false, drop =1:6)
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
