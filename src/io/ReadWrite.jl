"""
    PopGen.read(infile::String; kwargs...)
Wraps the individual file importers to read a file in as a `PopData` object. File type is
inferred from the file extension (case insensitive): \n

| File Format         | Extensions             | Docstring     |
| :------------------ | :--------------------- | :------------ |
| delimited           | `.csv`, `.txt`, `.tsv` | `?delimited`  |
| genepop             | `.gen`, `.genepop`     | `?genepop`    |
| structure           | `.str`, `.structure`   | `?structure`  |
| plink               | `.bed`, `.ped`  | `?plink`  |
| variant call format (vcf) | `.vcf`, `.vcf.gz`| `?vcf`  |
| variant call format (bcf) | `.bcf`, `.bcf.gz`| `?bcf`  |

This function uses the same keyword arguments (and defaults) as the file importing
functions it wraps; please see their respective docstrings in the Julia help console.
for specific usage details (e.g. `?genepop`).


## Examples
```
PopGen.read("cavernous_assfish.gen", digits = 3)

PopGen.read("juglans_nigra.vcf")
```
"""
function read(infile::String; kwargs...)::PopData
    ext = splitext(infile)[2] |> lowercase
    
    if occursin(r".vcf$", infile) | occursin(r".vcf.gz$", infile)
        return vcf(infile; kwargs...)
    
    elseif occursin(r".bcf$", infile) | occursin(r".bcf.gz$", infile)
        return bcf(infile; kwargs...)
    
    elseif ext in [".gen", ".genepop"]
        return genepop(infile;kwargs...)

    elseif ext in [".csv", ".txt", ".tsv"]
        return delimited(infile; kwargs...)

    elseif ext in [".str", ".structure"]
        return structure(infile; kwargs...)

    elseif ext in [".bed", ".ped"]
       return plink(infile; kwargs...)

    else
        @error "File type not recognized by filename extension \n delimited: .csv | .tsv | .txt \n genepop: .gen | .genepop \n structure: .str | .structure \n variant call formant: .bcf | .vcf \n plink: .bed | .ped | .fam | .bim | .map"
    end
end

"""
    PopGen.write(data::PopData, filename::String, kwargs...)
    PopGen.write(data::PopData; filename::String, kwargs...)

Writes `PopData` to a specified file type inferred from the extension of `filename = ` (case insensitive). Additional keyword
arguments `kwargs...` are specific to the intended file type, and are listed in the docstrings of the specific
file writer with the format `?filetype`. For example, to find the appropriate keywords for a conversion
to Genepop format, call up the `?genepop` docstring.

| File Format | Extensions             | Docstring          |
| :---------- | :--------------------- | :----------------- |
| genepop     | `.gen`, `.genepop`     | ?genepop   |
| delimited   | `.csv`, `.txt`, `.tsv` | ?delimited |
| structure   | `.str`, `.structure`   | ?structure |
| plink       | `.ped`                 | ?plink     |
| baypass     | `.baypass`             | ?baypass   |

### Example
```
cats = @nancycats;
fewer_cats = omit(cats, name = samplenames(cats)[1:10]);
PopGen.write(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")
```
"""
function Base.write(data::PopData; filename::String, kwargs...)
    ext = split(filename, ".")[end] |> lowercase
    if ext in ["gen", "genepop"]
        genepop(data, filename = filename; kwargs...)
    elseif ext in ["str", "structure"]
        structure(data, filename = filename; kwargs...)
    elseif ext in ["csv", "txt", "tsv"]
        delimited(data, filename = filename; kwargs...)
    elseif ext == "baypass"
        baypass(data, filename = filename)
    else
        throw(ArgumentError("File type not recognized by filename extension. Please see the docstring"))
    end
end

function Base.write(data::PopData, filename::String, kwargs...)
    ext = split(filename, ".")[end] |> lowercase
    if ext in ["gen", "genepop"]
        genepop(data, filename = filename; kwargs...)
    elseif ext in ["str", "structure"]
        structure(data, filename = filename; kwargs...)
    elseif ext in ["csv", "txt", "tsv"]
        delimited(data, filename = filename; kwargs...)
    elseif ext == "baypass"
        baypass(data, filename = filename)
    else
        throw(ArgumentError("File type not recognized by filename extension. Please see the docstring"))
    end
end