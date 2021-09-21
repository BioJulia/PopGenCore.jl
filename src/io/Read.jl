"""
    PopGen.read(infile::String; kwargs...)
Wraps the individual file importers to read a file in as a `PopData` object. File type is
inferred from the file extension (case insensitive): \n

| File Format         | Extensions             | Docstring     |
| :------------------ | :--------------------- | :------------ |
| delimited           | `.csv`, `.txt`, `.tsv` | `?PopGen.delimited`  |
| genepop             | `.gen`, `.genepop`     | `?PopGen.genepop`    |
| structure           | `.str`, `.structure`   | `?PopGen.structure`  |
| plink               | `.bed`, `.ped`  | `?PopGen.plink`  |
| variant call format (vcf) | `.vcf`, `.vcf.gz`| `?PopGen.vcf`  |
| variant call format (bcf) | `.bcf`, `.bcf.gz`| `?PopGen.bcf`  |

This function uses the same keyword arguments (and defaults) as the file importing
functions it wraps; please see their respective docstrings in the Julia help console.
for specific usage details (e.g. `?PopGen.genepop`). Replace `PopGen` with `PopGenCore` if 
using `PopGenCore.jl` directly.


## Examples
```
PopGen.read("cavernous_assfish.gen", digits = 3)

PopGen.read("juglans_nigra.vcf")
```
"""
function read(infile::String; kwargs...)
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
    PopGen.write(data::PopData; filename::String, kwargs...)
Writes `PopData` to a specified file type inferred from the extension of `filename = ` (case insensitive). Additional keyword
arguments `kwargs...` are specific to the intended file type, and are listed in the docstrings of the specific
file writer with the format `?filetype`. For example, to find the appropriate keywords for a conversion
to Genepop format, call up the `?genepop` docstring. Replace `PopGen` with `PopGenCore` if 
using `PopGenCore.jl` directly.

| File Format | Extensions             | Docstring          |
| :---------- | :--------------------- | :----------------- |
| genepop     | `.gen`, `.genepop`     | ?PopGen.genepop   |
| delimited   | `.csv`, `.txt`, `.tsv` | ?PopGen.delimited |
| structure   | `.str`, `.structure`   | ?PopGen.structure |

### Example
```
cats = @nancycats;
fewer_cats = omit(cats, name = samples(cats)[1:10]);
PopGen.write(fewer_cats, filename = "filtered_nancycats.gen", digits = 3, format = "h")
```
"""
function write_to(data::PopData; filename::String, kwargs...)
    ext = split(filename, ".")[end] |> lowercase
    if ext in ["gen", "genepop"]
        genepop(data, filename = filename; kwargs...)
    elseif ext in ["str", "structure"]
        structure(data, filename = filename; kwargs...)
    elseif ext in ["csv", "txt", "tsv"]
        delimited(data, filename = filename; kwargs...)
    else
        @error "File type not recognized by filename extension. Please see the docstring"
    end
end