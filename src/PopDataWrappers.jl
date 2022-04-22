"""
    info(::PopData)
Display the metadata (PopDataInfo) of a PopData object.
"""
metadata(data::PopData) = data.metadata
info(data::PopData) = data.metadata

"""
    sampleinfo(::PopData)
Show the sample information found within the metadata of a `PopData` object. Returns a view of a dataframe
"""
sampleinfo(data::PopData) = data.sampleinfo


"""
    genodata(::PopData)
Method to show the genotype information of a `PopData` object. Returns a view of a dataframe. 
"""
genodata(data::PopData) = data.genodata


"""
    locusinfo(::PopData)
Show the locus information found within the metadata of a `PopData` object. Returns a view of a dataframe
"""
locusinfo(data::PopData) = data.locusinfo

"""
    locationdata(data::PopData)
View the longitude and latitude data in a `PopData` object. Returns a table
derived from the PopData. Changes made to this table will not alter the source
`PopData` object.

Use `locations!` to add spatial data to a `PopData` object.
"""
function locationdata(data::PopData)
    if :longitude ∉ propertynames(data.sampleinfo) && :latitude ∉ propertynames(data.sampleinfo) 
        throw(ArgumentError(":longitude and :latitude columns not present in metadata."))
    else
        @view data.sampleinfo[!, [:longitude, :latitude]]
    end
end

"""
    loci(data::PopData)
Returns an array of strings of the loci names in a `PopData` object.
"""
function loci(data::PopData)
    data.locusinfo.locus
end

"""
    samplenames(data::PopData)
View individual/sample names in a `PopData`
"""
function samplenames(data::PopData)
    #string.(unique(data.genodata.name))
    string.(data.sampleinfo.name)
end

"""
    genotypes(data::PopData, samplelocus::String)
Return a vector of all the genotypes of a sample (or locus) in a `PopData` object.
```
cats = @nancycats
genotypes(cats, "N115")
genotypes(cats, "fca8")

```
"""
function genotypes(data::PopData, samplelocus::AbstractString)
    if samplelocus ∈ samplenames(data)
        @view data.genodata[data.genodata.name .== samplelocus, :genotype]
    elseif samplelocus ∈ loci(data)
        @view data.genodata[data.genodata.locus .== samplelocus, :genotype]
    else
        throw(ArgumentError("$samplelocus is not found in either samples or loci."))
    end
end

"""
    genotypes(data::PopData, samplelocus::Pair{String, String}) ::DataFrame
    genotypes(data::PopData, samplelocus::Pair{Vector{String}, String}) ::DataFrame
    genotypes(data::PopData, samplelocus::Pair{String, Vector{String}}) ::DataFrame
    
Return a genotype or dataframe of genotypes for one or more samples/loci 
in a `PopData` object. Uses the `Pair` notation of `samples => loci`.
### Examples
```
cats = @nancycats;
genotypes(cats, "N115" => "fca8")
genotypes(cats, ["N115", "N7"] => "fca8")
genotypes(cats, "N115" => ["fca8", "fca37"])
genotypes(cats, ["N1", "N2"] => ["fca8", "fca37"])
```
"""
function genotypes(data::PopData, samplocus::Pair{T,U})  where T<:AbstractString where U<:AbstractString
    sample, locus = samplocus
    @views data.genodata[(data.genodata.name .== sample) .& (data.genodata.locus .== locus), :genotype][1]
end

function genotypes(data::PopData, samplelocus::Pair{U, Vector{T}}) where T<:AbstractString where U<:AbstractString
    sample,locus = samplelocus
    @view data.genodata[(data.genodata.name .== sample) .& (data.genodata.locus .∈ Ref(locus)), :] 
end

function genotypes(data::PopData, samplelocus::Pair{Vector{U}, Vector{T}}) where T<:AbstractString where U<:AbstractString
    sample,locus = samplelocus
    @view data.genodata[(data.genodata.name .∈ Ref(sample)) .& (data.genodata.locus .∈ Ref(locus)), :] 
end

function genotypes(data::PopData, samplelocus::Pair{Vector{U}, T}) where T<:AbstractString where U<:AbstractString
    sample,locus = samplelocus
    @view data.genodata[(data.genodata.name .∈ Ref(sample)) .& (data.genodata.locus .== locus), :] 
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
