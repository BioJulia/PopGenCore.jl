export add_meta!, locations, locations!
export genotypes, get_genotypes, get_genotype
export populations, population, populations!, population!
export exclude, remove, omit, exclude!, remove!, omit!, keep, keep!, filter, filter!

# TODO make adding metata flexible to do inner joins


"""
    sampleinfo!(::PopData, metadata::Pair{Symbol, Vector}; categorical::Bool = false)
    sampleinfo!(::PopData, metadata::Pair{String, Vector}; categorical::Bool = false)
Add an additional sample information to `PopData` metadata. Mutates `PopData` in place. Metadata 
must be in the same order as the samples in `PopData.sampleinfo`.

#### Arguments
- `metadata` : A Pair of :ColumnName => [Values]

#### Keyword Arguments
- `categorical` : Boolean of whether the metadata being added is categorical aka "factors" (default: `false`)

## Example
```julia
cats = @nancycats
sampleinfo!(cats, :whiskerlength => rand(cats.metadata.samples))
sampleinfo!(cats, "tailcolor" => rand(["orange", "brown"], cats.metadata.samples), categorical = true)
cats
PopData{Diploid, 9 Microsatellite loci}
  Samples: 237
  Populations: 17
  Other Info: ["whiskerlength", "tailcolor"]
```
"""
function sampleinfo!(data::PopData, metadata::Pair{Symbol, T}; categorical::Bool = false) where T <: AbstractVector
    length(metadata[2]) != data.metadata.samples && error("Provided metadata vector (n = $(length(metadata[2]))) and samples in PopData (n = $(data.metadata.samples)) have different lengths")
    infotext = "\nAdding :$(metadata[1]) column to metadata.sampleinfo" 
    @info infotext
    if categorical == true
        insertcols!(data.metadata.sampleinfo, metadata[1] => PooledArray(metadata[2], compress = true))
    else
        insertcols!(data.metadata.sampleinfo, metadata[1] => metadata[2])
    end
    return
end

function sampleinfo!(data::PopData, metadata::Pair{String, T}; categorical::Bool = false) where T <: AbstractVector 
    sampleinfo!(data, Symbol(metadata[1]) => metadata[2], categorical = categorical)
end

"""
    locusinfo!(::PopData, metadata::Pair{Symbol, Vector}; categorical::Bool = false)
    locusinfo!(::PopData, metadata::Pair{String, Vector}; categorical::Bool = false)
Add an additional locus information to `PopData` metadata. Mutates `PopData` in place. Metadata 
must be in the same order as the samples in `PopData.locusinfo`.

#### Arguments
- `metadata` : A Pair of :ColumnName => [Values]

#### Keyword Arguments
- `categorical` : Boolean of whether the metadata being added is categorical aka "factors" (default: `false`)

## Example
```julia
cats = @nancycats
locusinfo!(cats, :quality => rand(cats.metadata.loci))
cats
PopData{Diploid, 9 Microsatellite loci}
  Samples: 237
  Populations: 17
  Other Info: ["quality"]
```
"""
function locusinfo!(data::PopData, metadata::Pair{Symbol, T}; categorical::Bool = false) where T <: AbstractVector
    length(metadata[2]) != data.metadata.loci && error("Provided metadata vector (n = $(length(metadata[2]))) and samples in PopData (n = $(data.metadata.loci)) have different lengths")
    infotext = "\nAdding :$(metadata[1]) column to metadata.locusinfo" 
    @info infotext
    if categorical == true
        insertcols!(data.metadata.locusinfo, metadata[1] => PooledArray(metadata[2], compress = true))
    else
        insertcols!(data.metadata.locusinfo, metadata[1] => metadata[2])
    end
    return
end

function locusinfo!(data::PopData, metadata::Pair{String, T}; categorical::Bool = false) where T <: AbstractVector 
    locusinfo!(data, Symbol(metadata[1]) => metadata[2], categorical = categorical)
end


"""
    locations(data::PopData)
View the longitude and latitude data in a `PopData` object. Returns a table
derived from the PopData. Changes made to this table will not alter the source
`PopData` object.

Use `locations!` to add spatial data to a `PopData` object.
"""
function locations(data::PopData)
    if :longitude ∉ propertynames(data.sampleinfo) && :latitude ∉ propertynames(data.sampleinfo) 
        throw(ArgumentError(":longitude and :latitude columns not present in metadata."))
    else
        @view data.sampleinfo[!, [:longitude, :latitude]]
    end
end


"""
    locations!(data::PopData; long::Vector{Float64}, lat::Vector{Float64})
Replaces existing `PopData` location data (longitude `long`, latitude `lat`).
Takes **decimal degrees** as a `Vector` of any `AbstractFloat`.
## Formatting requirements
- Decimal Degrees format: `-11.431`
- **Must** use negative sign `-` instead of cardinal directions
- Location data must be in the order that samples appear in your `PopData`
- Missing data should be coded as `missing` values of type `Missing` (can be accomplished with `replace!()`)
### Example
```
ncats = @nancycats ;
x = rand(237) ; y = rand(237)
locations!(ncats, long = x, lat = y)
```
"""
function locations!(data::PopData, long::Vector{Union{Missing,T}}, lat::Vector{Union{Missing,T}}) where T <: AbstractFloat
    long_len = length(long)
    lat_len = length(lat)
    long_len != lat_len && error("latitude ($lat_len) and longitude ($long_len) arrays not equal in length")
    long_len != length(data.sampleinfo.name) && error("lat/long array length ($long_len) and number of samples in PopData ($long_len) are not equal")

    data.sampleinfo.longitude = long
    data.sampleinfo.latitude = lat
    return
end

function locations!(data::PopData, long::Vector{T}, lat::Vector{T}) where T <: AbstractFloat
    # convert to the right type and use locations!()
    lat_adjust = lat |> Vector{Union{Missing, Float32}}
    long_adjust = long |> Vector{Union{Missing, Float32}}
    locations!(data, long = long_adjust, lat = lat_adjust)
end


"""
    locations!(data::PopData; long::Vector{String}, lat::Vector{String})
Replaces existing `PopData` location data (longitude `long`, latitude `lat`). Takes
**decimal minutes** or **degrees minutes seconds** format as a `Vector` of `String`. Recommended to use `CSV.read`
from `CSV.jl` to import your spatial coordinates from a text file.
## Formatting requirements
- Coordinates as a `String` separated by spaces (`"11 43 41"`) or colons (`"11:43:41"`)
- Must use negative sign (`"-11 43.52"`) or single-letter cardinal direction (`"11 43.52W"`)
- Missing data should be coded as the string `"missing"` (can be accomplished with `replace!()`)
- Can mix colons and spaces (although it's bad practice)
### NOTE
If you read in the coordinate data as 4 vectors (longitude degrees, longitude minutes, latitude degrees, latitude minutes),
then the easiest course of action would be to merge them into two vectors of strings
(one for longitude, one for latitude):
```
long_string = string.(lat_deg, " ", lat_min)
lat_string = string.(long_deg, " ", long_min)
```
and use these as inputs into `locations!`

### Example
```
ncats = @nancycats;
x = fill("11 22.33W", 237) ; y = fill("-41 31.52", 237)
locations!(ncats, long = x, lat = y)
```
"""
function locations!(data::PopData, long::Vector{String}, lat::Vector{String})
    long_len = length(long)
    lat_len = length(lat)
    lat_len != long_len && error("latitude ($lat_len) and longitude ($long_len) arrays not equal in length")
    lat_len != length(data.sampleinfo.name) && error("lat/long array length ($lat_len) and number of samples in PopData ($long_len) are not equal")
    println("Converting decimal minutes to decimal degrees")
    # convert coordinates to decimal degrees
    data.sampleinfo.longitude = convert_coord.(long)
    data.sampleinfo.latitude = convert_coord.(lat)
    return
end

function locations!(data::PopData; kwargs...)
    kwargs = Dict(kwargs)
    # check for matching lat and long keywords
    if all([haskey(kwargs, :lat), haskey(kwargs, :long)])
        locations!(data, kwargs[:long], kwargs[:lat])
    else
        error("keyword arguments \"lat\" and \"long\" must be supplied together")
    end
end


"""
    get_genotype(data::PopData; sample::String, locus::String)
Return the genotype of one sample at one locus in a `PopData` object.
### Example
```
cats = @nancycats;
get_genotype(cats, sample = "N115", locus = "fca8")
```
"""
function get_genotype(data::PopData; sample::String, locus::String)
    @views data.genodata[(data.genodata.name .== sample) .& (data.genodata.locus .== locus), :genotype][1]
end


"""
    get_genotypes(data::PopData, sample::String)
Return a vector of all the genotypes of a sample in a `PopData` object. To return a
single genotype at a locus, see `get_genotype`.
```
cats = @nancycats
get_genotypes(cats, "N115")
```
"""
function get_genotypes(data::PopData, sample::String)
    data.genodata[data.genodata.name .== sample, :genotype]
end

"""
    get_genotypes(::PopData; sample::Union{T, Vector{T}}, locus::Union{T, Vector{T}}) where T<:AbstractString
Return a table of the genotype(s) of one or more `samples` for one or more
specific `loci` (both as keywords) in a `PopData` object.
### Examples
```
cats = @nancycats;
get_genotypes(cats, name = "N115" , locus = "fca8")
get_genotypes(cats, name = ["N115", "N7"] , locus = "fca8")
get_genotypes(cats, name = "N115" , locus = ["fca8", "fca37"])
get_genotypes(cats, name = ["N1", "N2"] , locus = ["fca8", "fca37"])
```
"""
function get_genotypes(data::PopData; sample::Union{T, Vector{T}}, locus::Union{U, Vector{U}}) where T<:AbstractString where U<:AbstractString
    sample = typeof(sample) <: AbstractString ? [sample] : sample
    locus = typeof(locus) <: AbstractString ? [locus] : locus
    @view data.genodata[(data.genodata.name .∈ Ref(sample)) .& (data.genodata.locus .∈ Ref(locus)), :] 
end


"""
    genotypes(data::PopData, locus::Union{String, Symbol})
Convenience wrapper to return a vector of all the genotypes of a single locus

### Example
```
genotypes(@gulfsharks, "contig_475")
```
"""
function genotypes(data::PopData, locus::String)
    @view data.genodata[data.genodata.locus .== locus, :genotype]
end


"""
    populations(data::PopData; counts::Bool = false)
View unique population ID's and/or their counts in `PopData`.

- `counts` returns a dataframe of samples per `population` instead (default = `false`)
"""
@inline function populations(data::PopData; counts::Bool = false)
    if all(ismissing.(data.sampleinfo.population)) == true
        println("no population data present in PopData")
        return
    end
    uniq_pops = unique(data.sampleinfo.population)
    if counts == false
        return uniq_pops
    else
        pops = countmap(data.sampleinfo.population)
        return DataFrame(:population => uniq_pops, :count => [pops[i] for i in uniq_pops])
    end
end


"""
```
# Replace by matching
populations!(data::PopData, rename::Dict)
populations!(data::PopData, rename::Vector{String})
populations!(data::PopData, samples::Vector{String}, populations::Vector{String})
```
Multiple methods to rename or reassign population names to a `PopData`.

## Rename using a Dictionary
Rename existing population ID's of `PopData` using a `Dict` of
`population_name => replacement`
\n**Example**
```
potatopops = Dict("1" => "Idaho", "2" => "Russet")
populations!(potatoes, potatopops)
```

## Rename using a Vector of Strings
`Vector` of new unique population names. If the number of new names is equal to the number of current unique population names,
the method will intelligently rename the existing populations. If the number of new population names is equal to the number of samples,
the method will instead assign new population names to every sample in the order with which they appear in `PopData.metadata.sampleinfo`.
\n**Example**
```
# rename [2] existing populations
potatopops = ["Idaho", "Russet"]
populations!(potatoes, potatopops)

# assign new names to all [44] samples
poptatopops = repeat(["Idaho", "Russet"], inner = 22) ;
populations!(potatoes, potatopops)
```
## Reassign using samples and new population assignments
Completely reassign populations for each individual. Takes two vectors of strings
as input: one of the sample names, and the other with their new corresponding
population name. This can be useful to change population names for only some individuals.

\n**Example**
```
populations!(potatoes, ["potato_1", "potato_2"], ["north", "south"])
```

"""
function populations!(data::PopData, rename::Dict)
    msg = ""
    @inbounds for key in keys(rename)
        if key ∉ unique(data.sampleinfo.population)
            msg *= " population not found: \"$key\"\n"
        else
            replace!(data.sampleinfo.population, key => rename[key])
            replace!(data.genodata.population.pool, key => rename[key])
            data.genodata.population = PooledArray(data.genodata.population, compress = true)
            data.sampleinfo.population = PooledArray(data.sampleinfo.population, compress = true)
        end
    end
    msg != "" && printstyled("Warnings:", color = :yellow) ; print("\n"*msg)
    data.metadata.populations = length(unique(data.sampleinfo.population))
    return
end

function populations!(data::PopData, rename::Vector{String})
    current_popnames = unique(data.sampleinfo.population)
    if length(current_popnames) == length(rename)
        # infer that you want to replace existing unique names
        println(" Renaming unique populations")
        rn_dict = Dict(zip(current_popnames, rename))
        populations!(data, rn_dict)
    elseif length(rename) == data.metadata.samples
        # infer that you want to replace the population names for each sample
        println(" Assigning new population names to all samples")
        populations!(data, samples(data), rename)
    else 
        length(rename) != length(current_popnames) && throw(DimensionMismatch("Number of replacement names ($(length(rename))) do not match the number of current population names ($(length(current_popnames))) or number of samples ($(data.metadata.samples))"))
    end
    return
end

function populations!(data::PopData, samples::AbstractVector{T}, populations::AbstractVector{U}) where T<:AbstractString where U<:AbstractString
    nsample = length(samples)
    npops = length(populations)
    nsample != npops && throw(DimensionMismatch("Number of provided samples ($(nsample)) does not match the number of provided populations ($(npops))"))
    namescheck = symdiff(samples,intersect(samples, data.sampleinfo.name))
    if !isempty(namescheck)
        throw(ArgumentError("$(length(namescheck)) samples not found: $(join(namescheck, ", "))"))
    end
    for (i,samp) in enumerate(samples)
        @views data.sampleinfo[data.sampleinfo.name .== samp, :population] .= populations[i]
        @views data.genodata[data.genodata.name .== samp, :population] .= populations[i]
    end
    # drop old levels
    data.genodata.population = PooledArray(data.genodata.population, compress = true)
    data.sampleinfo.population = PooledArray(data.sampleinfo.population, compress = true)
    data.metadata.populations = length(unique(data.sampleinfo.population))
    return
end

##### Exclusion #####
"""
    exclude!(data::PopData, kwargs...)
Edit a `PopData` object in-place by excluding all occurences of the specified information.
The keywords can be used in any combination. Synonymous with `omit!` and `remove!`. All
values are converted to `String` for filtering, so `Symbol` and numbers will also work.
This can be considered a simpler and more rudimentary syntax for subsetting 
or filtering PopData.

### Keyword Arguments
#### `locus`
A `String` or `Vector{String}` of loci you want to remove from the `PopData`.
The keyword `loci` also works.

#### `population`
A `String` or `Vector{String}` of populations you want to remove from the `PopData`
The keyword `populations` also works.

#### `name`
A `String` or `Vector{String}` of samples you want to remove from the `PopData`
The keywords `names`, `sample`, and `samples` also work.

**Examples**
```
cats = @nancycats;
exclude!(cats, name = "N100", population = 1:5)
exclude!(cats, name = ["N100", "N102", "N211"], locus = ["fca8", "fca23"])
exclude!(cats, name = "N102", locus = :fca8, population = "3")
"""
function exclude!(data::PopData; population::Any = nothing, locus::Any = nothing, name::Any = nothing)
    filter_by = Dict{Symbol,Vector{String}}()

    if !isnothing(population)
        filter_by[:population] = typeof(population) <: AbstractArray ? string.(population) : [string(population)]
        err = filter_by[:population][filter_by[:population] .∉ Ref(unique(data.sampleinfo.population))]
        if length(err) > 0
            printstyled("Populations not found: ", bold = true)
            print("\"" * err[1] * "\"")
            if length(err) > 1
                [print(", \"$i\"") for i in err[2:end]] ; print("\n")
            end
            println()
        end
    end
    if !isnothing(locus)
        filter_by[:locus] = typeof(locus) <: AbstractArray ? string.(locus) : [string(locus)]
        err = filter_by[:locus][filter_by[:locus] .∉ Ref(loci(data))]
        if length(err) > 0
            printstyled("Loci not found: ", bold = true)
            print("\"" * err[1] * "\"")
            if length(err) > 1
                [print(", \"$i\"") for i in err[2:end]] ; print("\n")
            end
            println()
        end
    end
    if !isnothing(name)
        filter_by[:name] = typeof(name) <: AbstractArray ? string.(name) : [string(name)]
        err = filter_by[:name][filter_by[:name] .∉ Ref(data.sampleinfo.name)]
        if length(err) > 0
            printstyled("Samples not found: ", bold = true)
            print("\"" * err[1] * "\"")
            if length(err) > 1
                [print(", \"$i\"") for i in err[2:end]] ; print("\n")
            end
            println()
        end
    end

    filter_keys = Symbol.(keys(filter_by))

    if length(filter_keys) == 1
        filter!(data, filter_keys[1] => x -> x ∉ filter_by[filter_keys[1]])
    elseif length(filter_keys) == 2
        filter!(data, [filter_keys[1], filter_keys[2]] => (x,y) -> x ∉ filter_by[filter_keys[1]] && y ∉ filter_by[filter_keys[2]])
    elseif length(filter_keys) == 3
        filter!(data, [filter_keys[1], filter_keys[2], filter_keys[3]] => (x,y,z) -> x ∉ filter_by[filter_keys[1]] && y ∉ filter_by[filter_keys[2]] && z ∉ filter_by[filter_keys[3]])
    else
        throw(ArgumentError("Please specify at least one filter parameter of population, locus, or name"))   
    end

    return
end

const omit! = exclude!
const remove! = exclude!

"""
    exclude(data::PopData, kwargs...)
Returns a new `PopData` object excluding all occurrences of the specified keywords.
The keywords can be used in any combination. Synonymous with `omit` and `remove`. All
values are converted to `String` for filtering, so `Symbol` and numbers will also work.
This can be considered a simpler and more rudimentary syntax for subsetting 
or filtering PopData.

### Keyword Arguments
#### `locus`
A `String` or `Vector{String}` of loci you want to remove from the `PopData`.

#### `population`
A `String` or `Vector{String}` of populations you want to remove from the `PopData`.

#### `name`
A `String` or `Vector{String}` of samples you want to remove from the `PopData`.

**Examples**
```
cats = @nancycats;
exclude(cats, name = "N100", population = 1:5)
exclude(cats, name = ["N100", "N102", "N211"], locus = ["fca8", "fca23"])
exclude(cats, name = "N102", locus = :fca8, population = "3")
```
"""
function exclude(data::PopData; population::Any = nothing, locus::Any = nothing, name::Any = nothing)
    tmp = copy(data)
    exclude!(tmp; population = population, locus = locus, name = name)
    return tmp
end

const omit = exclude
const remove = exclude

"""
    keep!(data::PopData, kwargs...)
Edit a `PopData` object in-place by keeping only the occurrences of the specified keywords.
If using multiple fields, they will be chained together as "`or`" statements.
All values are converted to `String` for filtering, so `Symbol` and numbers will also work.
This can be considered a simpler and more rudimentary syntax for subsetting 
or filtering PopData.

### Keyword Arguments
#### `locus`
A `String` or `Vector{String}` of loci you want to keep in the `PopData`.

#### `population`
A `String` or `Vector{String}` of populations you want to keep in the `PopData`.

#### `name`
A `String` or `Vector{String}` of samples you want to keep in the `PopData`.

**Examples**

```
cats = @nancycats;
keep!(cats, population = 1:5)

# keep 4 populations and 3 specific samples
keep!(cats, name = ["N100", "N102", "N211"])

# keep 2 loci, 2 populations, and 10 specific individuals
keep!(cats, locus = [:fca8, "fca37"], population = [7,8], name = samples(cats)[1:10])
```
"""
function keep!(data::PopData; population::Any = nothing, locus::Any = nothing, name::Any = nothing)
    filter_by = Dict{Symbol,Vector{String}}()

    if !isnothing(population)
        filter_by[:population] = typeof(population) <: AbstractArray ? string.(population) : [string(population)]
        err = filter_by[:population][filter_by[:population] .∉ Ref(unique(data.sampleinfo.population))]
        if length(err) > 0
            printstyled("Populations not found: ", bold = true)
            print("\"" * err[1] * "\"")
            if length(err) > 1
                [print(", \"$i\"") for i in err[2:end]] ; print("\n")
            end
            println()
        end
    end
    if !isnothing(locus)
        filter_by[:locus] = typeof(locus) <: AbstractArray ? string.(locus) : [string(locus)]
        err = filter_by[:locus][filter_by[:locus] .∉ Ref(loci(data))]
        if length(err) > 0
            printstyled("Loci not found: ", bold = true)
            print("\"" * err[1] * "\"")
            if length(err) > 1
                [print(", \"$i\"") for i in err[2:end]] ; print("\n")
            end
            println()
        end
    end
    if !isnothing(name)
        filter_by[:name] = typeof(name) <: AbstractArray ? string.(name) : [string(name)]
        err = filter_by[:name][filter_by[:name] .∉ Ref(data.sampleinfo.name)]
        if length(err) > 0
            printstyled("Samples not found: ", bold = true)
            print("\"" * err[1] * "\"")
            if length(err) > 1
                [print(", \"$i\"") for i in err[2:end]] ; print("\n")
            end
            println()
        end
    end

    filter_keys = Symbol.(keys(filter_by))

    if length(filter_keys) == 1
        filter!(data, filter_keys[1] => x -> x ∈ filter_by[filter_keys[1]])
    elseif length(filter_keys) == 2
        filter!(data, [filter_keys[1], filter_keys[2]] => (x,y) -> x ∈ filter_by[filter_keys[1]] || y ∈ filter_by[filter_keys[2]])
    elseif length(filter_keys) == 3
        filter!(data, [filter_keys[1], filter_keys[2], filter_keys[3]] => (x,y,z) -> x ∈ filter_by[filter_keys[1]] || y ∈ filter_by[filter_keys[2]] || z ∈ filter_by[filter_keys[3]])
    else
        throw(ArgumentError("Please specify at least one filter parameter of population, locus, or name"))   
    end

    return
end


"""
    keep(data::PopData, kwargs...)
Returns a new `PopData` object keeping only the occurrences of the specified keyword.
Unlike `exclude()`. only one keyword can be used at a time. All values are 
converted to `String` for filtering, so `Symbol` and numbers will also work.
This can be considered a simpler and more rudimentary syntax for subsetting 
or filtering PopData.

### Keyword Arguments
#### `locus`
A `String` or `Vector{String}` of loci you want to keep in the `PopData`.

#### `population`
A `String` or `Vector{String}` of populations you want to keep in the `PopData`.

#### `name`
A `String` or `Vector{String}` of samples you want to keep in the `PopData`.

**Examples**
```
cats = @nancycats;
keep(cats, population = 1:5)
# equivalent to cats[cats.genodata.population .∈ Ref(string.(1:5)), :]

keep(cats, name = ["N100", "N102", "N211"])
# equivalent to cats[cats.genodata.name .∈ Ref(["N100", "N102", "N211"]), :]

keep(cats, locus = [:fca8, "fca37"])
# equivalent to cats[cats.genodata.locus .∈ Ref(["fca8", "fca37"]), :]
```
"""
function keep(data::PopData; population::Any = nothing, locus::Any = nothing, name::Any = nothing)
    tmp = copy(data)
    keep!(tmp; population = population, locus = locus, name = name)
    return tmp
end

"""
    filter(data::PopData, args...)
A drop-in replacement for DataFrames.filter where `PopData` is the first
argument and the filtering conditions are the second argument. Returns a
new `PopData`. **Note** the argument order is opposite of that from DataFrames.jl.

**Example**
```
x = @nancycats ;

y = filter(x, :name => i -> i ∈ samples(x)[1:10]) ;

show(x)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 237
  Populations: 17

show(y)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 10
  Populations: 1
```
"""
function filter(data::PopData, args...)
    out = copy(data) 
    filter!(out, args...)
    return out
end

"""
    filter(data::PopData, args...)
A drop-in replacement for the DataFrames.filter! where `PopData` is the first
argument and the filtering conditions are the second argument. Mutates the 
`PopData` in place and returns it. **Note** the argument order is opposite of that from DataFrames.jl.

**Example**
```
x = @nancycats ;

filter!(x, :name => i -> i ∈ samples(x)[1:10]) ;

show(x)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 10
  Populations: 1
```
"""
function Base.filter!(data::PopData, args...)
    geno = filter!(args..., data.genodata)
    transform!(
        geno,
        1 => (i -> PooledArray(i, compress = true)) => :name,
        2 => (i -> PooledArray(i, compress = true)) => :population,
        3 => (i -> PooledArray(i, compress = true)) => :locus,
        4
    )
    #=
    if intersect(data.sampleinfo.name, geno.name.pool) != data.sampleinfo.name
        filter!(:name => x -> x ∈ geno.name.pool, data.metadata)
    end
    =#
    PopDataInfo!(data)
    return data
end