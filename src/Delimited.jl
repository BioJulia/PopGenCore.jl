export delimited
"""
    delimited(infile::String; delim::Union{Char,String,Regex} = "auto", digits::Int64 = 3, silent::Bool = false)
Load a delimited-type file into memory as a PopData object.

### Arguments
- `infile` : path to file

### Keyword Arguments
- `delim` : delimiter characters. By default uses auto-parsing of `CSV.File`
- `digits` : number of digits denoting each allele (default: `3`)
- `silent`   : whether to print file information during import (default: `true`)
- `allow_monomorphic::Bool` : whether to keep monomorphic loci in the dataset (default: `false`)

### File formatting:
- The first row is column names (names don't matter)
- The columns must be in this order:
    1. sample name
    2. population name
    3. longitude
    4. latitude
    5. locus_1 genotypes
    6. locus_2 genotypes
    7. etc...

#### Missing data
##### Genotypes
Missing genotypes can be formatted as all-zeros `000000`, left empty, or negative-nine `-9`

##### Location data
If location data is missing for a sample (which is ok!), make sure the value is
blank, otherwise there will be transcription errors! (look at line 3 in the example below)

## Example
```
lizardsCA = delimited("CA_lizards.csv", digits = 3);
```

### Formatting example:
```
name,population,long,lat,Locus1,Locus2,Locus3 
sierra_01,mountain,11.11,-22.22,001001,-9,001001
sierra_02,mountain,11.12,-22.21,001001,001001,001002
snbarb_01,coast,,,001001,000000,001002
snbarb_02,coast,11.14,-22.24,001001,001001,001001
snbarb_03,coast,11.15,,001002,001001,001001
```
"""
function delimited(
    infile::String;
    delim::Union{Char,String,Regex} = "auto",
    digits::Int = 3,
    silent::Bool = false,
    allow_monomorphic::Bool = false
    )

    dlm = delim == "auto" ? nothing : delim
    file_parse = CSV.File(infile, delim = dlm, missingstrings = ["-9", ""]) |> DataFrame
    locinames = names(file_parse)[5:end]
    meta = select(file_parse, 1:4)
    # force strings for this field
    meta.population = string.(meta.population)
    rename!(meta, [:name, :population, :longitude, :latitude])
    select!(file_parse, Not(3:4))
    geno_parse = DataFrames.stack(file_parse, DataFrames.Not(1:2))
    rename!(geno_parse, [:name, :population, :locus, :genotype])
    
    select!(geno_parse,
        :name => (i -> PooledArray(i, compress = true)) => :name,
        :population => (i -> PooledArray(string.(i), compress = true)) => :population,
        :locus => (i -> PooledArray(i, compress = true)) => :locus, 
        :genotype
    )
    
    try
        geno_parse.genotype = phase.(geno_parse.genotype, Int8, digits)
    catch
        geno_parse.genotype = phase.(geno_parse.genotype, Int16, digits)
    end
    
    if !silent
        @info "\n $(abspath(infile))\n data: loci = $(length(locinames)), samples = $(length(meta[!, 1])), populations = $(length(unique(meta[!,2])))"
        println()
    end

    ploidy = DataFrames.combine(
        groupby(dropmissing(geno_parse), :name),
        :genotype => find_ploidy => :ploidy
    ).ploidy
    DataFrames.insertcols!(meta, 3, :ploidy => ploidy)
    pd_out = PopData(meta, geno_parse)
    !allow_monomorphic && drop_monomorphic!(pd_out, silent = silent)
    return pd_out
end


"""
    delimited(data::PopData; filename::String, delim::String = ",", digits::Integer = 3, format::String = "wide", miss::Int = 0)
Write PopData to a text-delimited file.
### Keyword Arguments
- `filename`: a `String` of the output filename
- `digits` : an `Integer` indicating how many digits to format each allele as (e.g. `(1, 2)` => `001002` for `digits = 3`)
- `format` : a `String` indicating whether to output in`"wide"` or `"long"` (aka `"tidy"`) format
  - `wide` : the standard format CSV for importing into PopGen.jl
  - `long` : the `loci` table with `longitude` and `latitude` columns added
- `delim` : the `String` delimiter to use for writing the file
- `miss` : an `Integer` for how you would like missing values written 
    - `0` : As a genotype represented as a number of zeroes equal to `digits × ploidy` like `000000` (default) 
    - `-9` : As a single value `-9`

### Example
```julia
cats = @nancycats;
fewer_cats = omit(cats, name = samples(cats)[1:10]);
delimited(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "wide", delim = " ")
```
"""
function delimited(data::PopData; filename::String, delim::String = ",", digits::Integer = 3, format::String = "wide", miss::Int = 0)
    # create empty String column
    unphased_df = insertcols!(copy(data.genodata), :string_geno => fill("", length(data.genodata.name)))
    # grouping to facilitate pulling out ploidy values for coding missing genotypes
    grp_df = groupby(unphased_df, :name)
    for sample in grp_df
        samplename = sample.name[1]
        sample_ploidy = convert(Int, data.metadata.ploidy[data.metadata.name .== samplename][1])
        sample.string_geno .= unphase.(sample.genotype, digits = digits, ploidy = sample_ploidy, miss = miss)
    end
    # pop out original genotype column
    select!(unphased_df, 1, 2, 3, 5)
    if format == "wide"
        wide_df = DataFrames.unstack(unphased_df, :locus, :string_geno)
        insertcols!(wide_df, 3, :longitude => data.metadata.longitude, :latitude => data.metadata.latitude)
        CSV.write(filename, wide_df, delim = delim) ;
        return
    else
        out_df = sort(
            transform(
                unphased_df, 1:2,
                3 => (i -> Vector{Union{Float32, Missing}}(undef, length(i))) => :longitude,
                4 => (i -> Vector{Union{Float32, Missing}}(undef, length(i))) => :latitude, 3:4
                )
                    )
        for i in 1:length(data.metadata.name)
            out_df[out_df.name .== data.metadata.name[i], :longitude] .= data.metadata.longitude[i]
            out_df[out_df.name .== data.metadata.name[i], :latitude] .= data.metadata.latitude[i]
        end
        CSV.write(filename, out_df, delim = delim) ;
        return
    end
end