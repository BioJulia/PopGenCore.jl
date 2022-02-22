"""
    genepop(infile::String; kwargs...)
Load a Genepop format file into memory as a PopData object.
- `infile::String` : path to Genepop file

### Keyword Arguments
- `digits::Integer`: number of digits denoting each allele (default: `3`)
- `popsep::String` : word that separates populations in `infile` (default: "POP")
- `silent::Bool`   : whether to print file information during import (default: `false`)
- `allow_monomorphic::Bool` : whether to keep monomorphic loci in the dataset (default: `false`)


### File must follow standard Genepop formatting:
- First line is a comment (and skipped)
- Loci are listed after first line as one-per-line without commas or in single comma-separated row
- A line with a particular keyword (default `POP`) must delimit populations
- Sample name is immediately followed by a *comma*
- File is *tab or space delimted* (but not both!)

### Genepop file example:
```
wasp_hive.gen: Wasp populations in New York
Locus1
Locus2
Locus3
pop
Oneida_01,  250230  564568  110100
Oneida_02,  252238  568558  100120
Oneida_03,  254230  564558  090100
pop
Newcomb_01, 254230  564558  080100
Newcomb_02, 000230  564558  090080
Newcomb_03, 254230  000000  090100
Newcomb_04, 254230  564000  090120
```
## Example
```
waspsNY = genepop("wasp_hive.gen", digits = 3, popsep = "pop")
```
"""
function genepop(
    infile::String;
    digits::Int = 3,
    popsep::String = "POP",
    silent::Bool = false,
    allow_monomorphic::Bool = false
)
    isfile(infile) || throw(ArgumentError("$infile not found."))
    # open the file as lines of strings to suss out loci names, pop idx, and popcounts
    gpop_readlines = readlines(infile)

    # find the #samples per population
    pop_idx = Vector{Int}()
    for (line, gtext) in enumerate(gpop_readlines)
        gtext == popsep && push!(pop_idx, line)
    end
    length(pop_idx) == 0 &&
    error("No populations found in $infile using separator \"$popsep\". Please check the spelling and try again.")

    # create a theoretical place were the last popsep would be (preserve last population for counting)
    push!(pop_idx, countlines(infile) + 1)
    popcounts = (pop_idx[2:end] .- pop_idx[1:end-1]) .- 1

    # check for the delimiter type in the first sample record
    firstrecord = gpop_readlines[pop_idx[1]+1]
    if occursin("\t", firstrecord) & occursin(" ", firstrecord)
        error("$infile contains both tab and space delimiters. Please format the file so it uses either one or the other.")
    elseif occursin("\t", firstrecord)
        delim = "\t"
        delim_txt = "tab"
    elseif occursin(" ", firstrecord)
        delim = " "
        delim_txt = "space"
    else
        error("Please format $infile to be either tab or space delimited")
    end

    # check for horizontal formatting, where popsep would appear on the third line
    if pop_idx[1] <= 3
        # second line should have all the loci
        locus_name_raw = replace(gpop_readlines[2], "." => "_")
        locinames = strip.(split(locus_name_raw |> join, ","))
        format = "horizontal"
    else
        #standard  vertical formatting
        locinames = gpop_readlines[2:pop_idx[1]-1]
        format = "vertical"
    end

    if !silent
        @info "\n $(truncatepath(abspath(infile)))\n formatting: delimiter = $(delim_txt), loci = $(format)\n data: loci = $(length(locinames)), samples = $(sum(popcounts)), populations = $(length(popcounts))"
        println()
    end

    # load in samples and genotypes
    pushfirst!(locinames, "name")

    geno_parse = CSV.read(
        infile,
        DataFrame,
        delim = delim,
        header = locinames,
        skipto = pop_idx[1] + 1,
        comment = popsep,
        missingstring = ["-9", ""],
        normalizenames = true,
        ignorerepeated = true
    )

    popnames = string.(1:length(popcounts))
    popnames = fill.(popnames,popcounts) |> Base.Iterators.flatten |> collect
    insertcols!(geno_parse, 2, :population => popnames)
    geno_parse.name .= strip.(geno_parse.name, ',')
    # wide to long format
    geno_parse = DataFrames.stack(geno_parse, DataFrames.Not([1,2]))
    rename!(geno_parse, [:name, :population, :locus, :genotype])
    select!(
        geno_parse, 
        :name => (i -> PooledArray(i, compress = true)) => :name, 
        :population => (i -> PooledArray(i, compress = true)) => :population, 
        :locus => (i -> PooledArray(Array(i), compress = true)) => :locus, 
        :genotype
    )

    # try to compress the alleles into Int8 (snps) or In16 (msat)
    try
        geno_parse.genotype = PooledArray(phase.(geno_parse.genotype, Int8, digits), compress = true)
    catch
        geno_parse.genotype = phase.(geno_parse.genotype, Int16, digits)
    end
    pd_out = PopData(geno_parse)
    !allow_monomorphic && dropmonomorphic!(pd_out, silent = silent)
    return pd_out
end

"""
    genepop(data::PopData; filename::String = "output.gen", digits::Int = 3, format::String = "vertical", miss::Int = 0)
Writes a `PopData` object to a Genepop-formatted file.
- `data`: the `PopData` object you wish to convert to a Genepop file
### keyword arguments
- `filename`: a `String` of the output filename
- `digits` : an `Integer` indicating how many digits to format each allele as (e.g. `(1, 2)` => `001002` for `digits = 3`)
- `format` : a `String` indicating whether loci should be formatted
    - vertically (`"v"` or `"vertical"`)
    - hortizontally (`"h"`, or `"horizontal"`)
    - Genepop Isolation-By-Distance (`"ibd"`) where each sample is a population with long/lat data prepended
- `miss` : an `Integer` for how you would like missing values written 
    - `0` : As a genotype represented as a number of zeroes equal to `digits × ploidy` like `000000` (default) 
    - `-9` : As a single value `-9`

## Example
```julia
cats = @nancycats;
fewer_cats = omit(cats, name = samplenames(cats)[1:10]);
genepop(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")
```
"""
function genepop(data::PopData; filename::String = "output.gen", digits::Int = 3, format::String = "vertical", miss::Int = 0)
    outfile = open(filename, "w") 
    println(outfile, "generated from PopData by PopGen.jl")
    if format in ["h", "horizontal"]
        join(outfile, loci(data), ",")
    else
        join(outfile, loci(data), "\n")
    end
    print(outfile, "\n")
    if lowercase(format) != "ibd"
        pops = Vector{String}()
        for (keys, sample) in pairs(groupby(data.genodata, [:name, :population]))
            if keys.population ∉ pops
                push!(pops, keys.population)
                println(outfile, "POP")
            end
            samplename = keys.name
            sample_ploidy = convert(Int, data.sampleinfo.ploidy[data.sampleinfo.name .== samplename][1])
            format_geno = unphase.(sample.genotype, digits = digits, ploidy = sample_ploidy, miss = miss)
            join(outfile, vcat(samplename * ",", format_geno), "\t")
            print(outfile, "\n")
        end
    else
        for (keys, sample) in pairs(groupby(data.genodata, :name))
            samplename = keys.name
            sample_ploidy = convert(Int, data.sampleinfo.ploidy[data.sampleinfo.name .== samplename][1])
            println(outfile, "POP")
            long = data.sampleinfo[data.sampleinfo.name .== keys.name, :longitude][1]
            lat = data.sampleinfo[data.sampleinfo.name .== keys.name, :latitude][1]
            coords = join(string.([long, lat]), "\t") * ","
            format_geno = unphase.(sample.genotype, digits = digits, ploidy = sample_ploidy, miss = miss)
            join(outfile, vcat(coords, format_geno), "\t")
            print(outfile, "\n")
        end
    end
    close(outfile)
end

precompile(genepop, (String,))
precompile(genepop, (PopData,))