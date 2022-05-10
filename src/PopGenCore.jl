module PopGenCore

using CSV, DataFrames, PooledArrays, StaticArrays
using CodecZlib, VariantCallFormat
using StatsBase: countmap, proportionmap
using Random: shuffle, shuffle!

include("PopData.jl")
export PopObj, PopData, PopDataInfo, show, Genotype, GenoArray, SNP, MSat
export getindex, getproperty

include("PopDataWrappers.jl")
export metadata, info, sampleinfo, genodata, locusinfo, samplenames
export locationdata, loci, samples, populations
export genotype, genotypes

## Utilities
include("Utils/GeneralUtils.jl")
export copy, size, sort, convertcoord
export dropmonomorphic, dropmonomorphic!
export dropmultiallelic, dropmultiallelic!

include("Utils/GenotypeUtils.jl")
export allelecount, alleles, uniquealleles
export locidataframe, locimatrix, phasedmatrix

include("Utils/ioUtils.jl")
export phase, unphase, findploidy

include("Utils/MathUtils.jl")
export reciprocal, reciprocalsum, countnonzeros

include("Utils/MissingUtils.jl")
export nonmissing, nonmissings

include("AlleleFreq.jl")
export allelefreq, allelefreq_vec, avg_allelefreq

include("GenoFreq.jl")
export genofreq, genofreq_expected, genocount_observed, genocount_expected

include("Conditionals.jl")
export ishom, ishet, _ishom, _ishet, isbiallelic, isbinary

include("Permutations.jl")
export permuteloci!, permutesamples!, permutegenotypes!, permutealleles!
export strictshuffle, strictshuffle!

include("Iterators.jl")
export skipinf, skipnan, skipinfnan, partitionarray, pairwisepairs, simpairs

include("Manipulate.jl")
export locationdata!, populations!, sampleinfo!, locusinfo!
export exclude, remove, omit, exclude!, remove!, omit!, keep, keep!, filter, filter!

## File IO
include("io/Delimited.jl")
export delimited

include("io/Genepop.jl")
export genepop

include("io/Structure.jl")
export structure

include("io/VariantCall.jl")
export bcf, vcf

include("io/ReadWrite.jl")
##
include("Datasets.jl")
export @nancycats, @gulfsharks, dataset

end