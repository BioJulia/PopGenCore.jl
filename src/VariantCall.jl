export bcf, vcf

"""
    bcf(infile::String; ; rename_loci::Bool, silent::Bool, allow_monomorphic::Bool)
Load a BCF file into memory as a PopData object. Population information needs to be provided separately. 
- `infile` : path to BCF file (can be gzipped)

**Keyword Arguments**
- `rename_loci` : true/false of whether to simplify loci names to "snp_#" (default: `false`)
- `allow_monomorphic` : true/false of whether to keep monomorphic loci (default: `false`)
- `silent`: true/false of whether to print extra file information (default: `false`).
Alleles are recoded according to the following schema:


| **Base**   |  A   |  T   |  C   |  G   |
| :--------  | :--: | :--: | :--: | :--: |
| **Allele** |  1   |  2   |  3   |  4   |


### Mixed-ploidy data
If importing mixed-ploidy data (such as poolseq), you will need to perform an additional
step to convert the genotype column into the correct `GenoArray` type:
```julia
julia> mydata = bcf("path/to/file.bcf", silent = true, rename_loci = true) ;

julia> mydata.genodata.genotype =  mydata.genodata.genotype |> Array{Union{Missing, NTuple}}
```
"""
bcf(infile::String; kwargs...) = error("Please load in VariantCallFormat.jl with \`using VariantCallFormat\`")

### VCF parsing ###

"""
    vcf(infile::String; ; rename_loci::Bool, silent::Bool, allow_monomorphic::Bool)
Load a VCF file into memory as a PopData object. Population information needs to be provided separately. 
- `infile` : path to VCF file (can be gzipped)

**Keyword Arguments**
- `rename_loci` : true/false of whether to simplify loci names to "snp_#" (default: `false`)
- `allow_monomorphic` : true/false of whether to keep monomorphic loci (default: `false`)
- `silent`: true/false of whether to print extra file information (default: `false`).
Alleles are recoded according to the following schema:


| **Base**   |  A   |  T   |  C   |  G   |
| :--------  | :--: | :--: | :--: | :--: |
| **Allele** |  1   |  2   |  3   |  4   |


### Mixed-ploidy data
If importing mixed-ploidy data (such as poolseq), you will need to perform an additional
step to convert the genotype column into the correct `GenoArray` type:
```julia
julia> mydata = vcf("path/to/file.vcf", silent = true, rename_loci = true) ;

julia> mydata.genodata.genotype =  mydata.genodata.genotype |> Array{Union{Missing, NTuple}}

```
"""
vcf(infile::String; kwargs...) = error("Please load in VariantCallFormat.jl with \`using VariantCallFormat\`")

