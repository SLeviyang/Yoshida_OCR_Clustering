# Workflow Steps

All data is stored in `s3://yoshida-ataqseq/`

## Step 1: Download fastq from SRA to S3

Folder : `fastq/`

## Step 2: Generate bam files by aligning fastq to mm10 using bowtie2

Folder : `bam/`

## Step 3: Call peaks for each bam file using macs2

FOLDER: `peaks/`

INPUT: 

1. bam files for a given cell type

OUTPUT:

1. peaks based on each replicate or on 2 pseudoreplicates if there is only 1 replicate for the cell type
2. peaks based on pooled replicates

## Step 4: Apply idr algorithm to find replicatable peaks

FOLDER: `idr/`

* Uses pooled peak calls as oracle list and then chooses the largest set of peaks called using idr at significance 0.05 over all pairs of peaks.

INPUT:

1. peak calls for replicates (possibly pseudo)
2. peak calls for pooled replicates

OUTPUT:

1. A bed+6 file genereted by the idr program



## Step 5: Find ISRE peaks within an idr peak

* access the cell type bed file from `idr/` folder in S3
* filter all isre to isre interesecting an idr peak 

FOLDER: `ISRE/`

INPUT: 

1. `idr/` files
2. bed file of all isre over mm10 (produced by the R script `make_all_ISRE_bed.R` in the step 5 folder)

OUTPUT: 

bed file with 8 columns (chr, chrStart, chrEnd, peakName, score, strand, seq, revC).  

## Step 6a: Use all ISRE peaks to define a master list of ISRE peaks found in at least one cell type

* download all the isre files from `ISRE/` folder in S3
* merges filtered isre across the cell types keeping only unique copies

INPUT:

1. a bed file of all isre in mm10
2. the `ISRE/` bed files  generated in the previous step for isre within each cell type

OUTPUT:

1. `master_isre.bed` file containing all isre within 10 base pair of a idr peak within at least one cell type
2. `master_isre_nonoverlapping.bed` file, which is same as `master_isre.bed` except that isre are filtered so that they don't overlap

## Step 6b: Use all idr peaks to define a master list of idr peaks found in at least one cell type

* download all the idr files from `idr/` folder in S3
* define a idr peak as the summit +/- 250 base pair
* merges filtered isre across the cell types and choose a greedy subset that of non-intersecting idr windows

INPUT:

1. the `idr/` bed files  

OUTPUT:

2. `master_idr_nonoverlapping.bed` file which contains columns (chr, chrStart, chrEnd, q)

## Step 7: Merge cell type ISRE to form a binary ISRE matrix

### Note: 

In previous workflow versions, I followed Yoshida and calculated read counts over master isre for each cell type.  The problem with this approach is the need to normalize both across cell types and across loci.

Observations

1. Quantile normalization of the cell types is problematic because there may not be a collection of loci that is implicitely the same across the cell types and acts as a baseline.
2. Normalization across loci has the issue of bias due to distance to TSS.

