export structure

"""
	phase_structure(datatype::DataType, args...)
Takes a DataType (such as `Int8`) and a series of integers to return
a sorted Tuple of those integers converted to that DataType. i.e. takes
a series of alleles and returns a genotype. Returns `missing` if args are
`missing`. Used internally in PopGen.structure file reader.
**Example**
```
phase_structure(Int8, 1,2,3,4,3,4,6,1)
(1, 1, 2, 3, 3, 4, 4, 6)
phase_structure(Int16, missing, missing)
missing
```
"""
function phase_structure(datatype::DataType, args...)
    all(ismissing.(args)) && return missing
    return Tuple(datatype.(sort([args...])))
end

"""
    structure(infile::String; kwargs...)
Load a Structure format file into memory as a PopData object.
- `infile::String` : path to Structure file

### Keyword Arguments
- `extracols::Integer`: how many additional optional columns there are beyond Stucture's POPDATA the reader needs to ignore (default: `0`)
    - these include POPFLAG, LOCDATA, or anything else you might have added
- `extrarows::Integer` : how many additional optional rows there are beyond the first row of locus names (default: `0`)
- `missingval::String`  : the value used to identify missing values in the data (default: `"-9"`)
- `silent::Bool`   : whether to print file information during import (default: `false`)
- `allow_monomorphic::Bool` : whether to keep monomorphic loci in the dataset (default: `false`)
- `faststructure::Bool`: whether the file is fastStructure format (default: `false`)

### File must follow this Structure format:
- the file is `tab` or `space` delimited **but not both**
- first row is locus names separated by the delimiter
    - leading/trailing whitespaces are tolerated
    - optional rows allowed **after** the locus names
- number of rows per sample = ploidy
    - e.g. if diploid, that sample would have 2 rows
    - multi-column variant not supported
    - all samples must have the same ploidy
- first data column is sample name
- second data column is population ID
    - optional columns allowed **after** the population ID (2nd) column
- remaining columns are the genotype for that individual for that locus

### Structure file example:
```
locus_1	locus_2	locus_3	locus_4	locus_5
walnut_01	1	-9	145	66	0	92
walnut_01	1	-9	-9	64	0	94
walnut_02	1	106	142	68	1	92
walnut_02	1	106	148	64	0	94
walnut_03	2	110	145	-9	0	92
walnut_03	2	110	148	66	1	-9
```

### fastStructure file format:
- the file is `tab` or `space` delimited **but not both**
- no first row of loci names
- number of rows per sample = ploidy
    - e.g. if diploid, that sample would have 2 rows
- first data column is sample name
- second data column is population ID
- remaining columns are the genotype for that individual for that locus
- usually, first 6 colums are empty (but not necessary)
- **no** extra rows or columns.
### fastStructure file example:
```
chestnut_01	1	-9	145	66	0	92
chestnut_01	1	-9	-9	64	0	94
chestnut_02	1	106	142	68	1	92
chestnut_02	1	106	148	64	0	94
chestnut_03	2	110	145	-9	0	92
chestnut_03	2	110	148	66	1	-9
```

## Example
```
walnuts = structure("juglans_nigra.str", extracols = 0, extrarows = 0)
```
"""
function structure(infile::String; silent::Bool = false, extracols::Int = 0, extrarows::Int = 0, allow_monomorphic::Bool = false, missingval::String = "-9", faststructure::Bool = false)
    # find the delimiter
    delimcheck = join(open(readlines, `head -n $(2+extrarows) $(infile)`))

    if occursin("\t", delimcheck) & occursin(" ", delimcheck)
        error("$infile contains both tab and space delimiters. Please format the file so it uses either one or the other.")
    elseif occursin("\t", delimcheck)
        delim = "\t"
    elseif occursin(" ", delimcheck)
        delim = " "
    else
        error("Please format $infile to be either tab or space delimited")
    end
    
    if faststructure == true
        data_row = 1
        first_row = split(strip(open(readline, infile)), delim)
        n_loci = length(first_row)
        locinames = ["locus_$i" for i in 1:n_loci-2]
    else
        data_row = 2 + extrarows
        locinames = string.(split(strip(open(readline, infile)), delim))
    end
    nloc = length(locinames)
    dropcols = extracols == 0 ? [1,2] : collect(1:2+extracols)
    names_parse = CSV.File(
        infile,
        delim = delim,
        header = false,
        datarow = data_row,
        missingstrings = [missingval],
        type = String,
        ignorerepeated = true,
        select = [1,2]
     ) |> DataFrame

    samplecount = countmap(names_parse[:,1])
    ploidy = unique(values(samplecount))
    # TODO make this flexible enough to take any ploidy
    length(ploidy) >1 && throw(error("Multiple ploidies detected among samples. Samples must all be of a single ploidy"))
    samplenames = unique(names_parse)  

    # read in the file as a table
    geno_parse = CSV.File(
        infile,
        delim = delim,
        header = false,
        datarow = data_row,
        missingstrings = [missingval],
        type = Int32,
        ignorerepeated = true,
        drop = dropcols
     ) |> Tables.matrix

    genos = Base.Iterators.partition(geno_parse, first(ploidy)) .|> sort!
    try
        genos = _SNP.(genos)
    catch
        genos = _MSat.(genos)
    end

    unique!(names_parse)
    rename!(names_parse, [:name, :population])
    nsamples = length(names_parse.name)
    
    #geno_parse.name .= replace.(geno_parse.name, "-" => "_")
    loci_df = DataFrame(
        :name => PooledArray(repeat(names_parse.name, nloc), compress = true),
        :population => PooledArray(repeat(names_parse.population, nloc), compress = true),
        :locus => PooledArray(repeat(locinames, inner = nsamples), compress = true),
        :genotype => genos,
    )

    # fix names, just in case
    
    
    #= create new dataframe to add phased genotypes to
    by_sample = groupby(geno_parse, :name)
    loci_df = DataFrame(:locus => locinames)
    for (key, eachsample) in pairs(by_sample)
        insertcols!(
            loci_df,    
            Symbol(key.name) => map(i -> phase_structure(markertype, i...), eachcol(eachsample)[3:end])
        )
    end

    # transpose the dataframe
    loci_df = permutedims(loci_df, 1, :name)
    =#

    # create the metadata from the original file info
    names_parse.ploidy = [samplecount[i] for i in names_parse.name]
    names_parse.longitude = Vector{Union{Missing, Float32}}(undef, nsamples) 
    names_parse.latitude = Vector{Union{Missing, Float32}}(undef, nsamples)

    if !silent
        @info "\n $(abspath(infile))\n data: loci = $(nloc), samples = $(nsamples), populations = $(length(unique(names_parse.population)))"
        println()
    end

    pd_out = PopData(names_parse, loci_df)
    !allow_monomorphic && drop_monomorphic!(pd_out) 

    return pd_out
end


"""
    structure(data::PopData; filename::String, faststructure::Bool, delim::String)
Write a `PopData` object to a Stucture format file
- `data`: the `PopData` object you wish to write to a Structure file
### keyword arguments
- `filename`: a `String` of the output filename
- `delim` : a `String` of either `"tab"` or `"space"` indicating the delimiter (default: `"tab"`)
- `faststructure`: true/false of whether the output should be formatted for fastStructure (default: `false`)

```
cats = @nancycats;
fewer_cats = omit(cats, name = samples(cats)[1:10]);
structure(fewer_cats, filename = "filtered_nancycats.str", faststructure = true)
```
"""
function structure(data::PopData; filename::String, faststructure::Bool = false, delim::String = "tab")
    # index both dataframes
    genos_gdf = groupby(data.genodata, :name)
    meta_gdf = groupby(data.metadatadata, :name)
    # get the sample names to iterate keys over
    idx = collect(samples(data))
    
    outfile = open(filename, "w")
    
    # check delimiter
    if delim == "tab"
        dlm = "\t"
    elseif delim == "space"
        dlm = " "
    else
        throw(ArgumentError("Please choose from either \"tab\" (default) or \"space\" delimiters."))
    end

    faststructure == false && println(outfile, join([i * dlm for i in loci(data)]))
    
    # remap populations as integers
    pops = unique(data.metadatadata.population)
    pop_mappings = Dict{String,Integer}()
    [pop_mappings[j] = i for (i,j) in enumerate(pops)]

    for sampl in idx
        ploid = meta_gdf[(name = sampl,)].ploidy |> first
        pop_id = meta_gdf[(name = sampl,)].population |> first
    
        # copy so as to not overwrite
        genos = copy(genos_gdf[(name = sampl,)].genotype)
    
        # replace missing with -9's
        miss_idx = findall(ismissing, genos)
        for i in miss_idx
            genos[i] = phase_structure(Int8, fill(-9, ploid)...)
        end

        # write the alleles to the file
        for allele in 1:ploid
            tmp_alleles = getindex.(genos, allele)
            allele_row = join([string(i, dlm) for i in tmp_alleles[1:end-1]])
            allele_row *= string(tmp_alleles[end])
            println(outfile, sampl, dlm, pop_mappings[pop_id], dlm, allele_row)
        end
    end
    close(outfile)
end
