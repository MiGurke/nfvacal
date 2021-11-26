# Nfvacal - nextflow variant calling pipeline using freebayes

 This pipeline calls variants from bam files of multiple individuals using [freebayes](https://github.com/freebayes/freebayes). It is adapted for usage on the MfN cluster.

Steps of the pipeline are:

1. Adding RG individual tags to the bam file to later be able to identify different individuals in the final cvf file.
2. Index bam files./
3. Splitting the scaffolds of the reference genome into sections of a given size.
4. Calling of variants for all individuals simultaneously, but in parallel for the sections of the scaffolds genome using freebayes.
5. Compressing the resulting vcf files.
6. Sorting the samples in the vcf file.
7. Indexing the files.
8. Merging of all chromosome vcf's into one final file.
9. Calculating statistics about the file using vcftools and vcflib.
10. Creating a html report from the stats. 



## Usage

```bash
nextflow run vacal.nf --bamdir /PATH/TO/BAMS/ --ref /PATH/TO/REFERENCE/file.fasta --outdir /PATH/TO/DIR/TO/PUT/OUTPUT/ --chr /PATH/TO/CHROMOSOME/file.list -with-report /PATH/TO/REPORT/file.html --split INT --poplist /PATH/TO/POP.LIST
```

**Parameters:**

* --bamdir : Path to directory in with sub folders for each individual. Each sub folder must hold one bam file.
* --ref : Path to the reference genome assembly fasta file.
* --outdir : Path to directory were the output of the pipeline will be written into.
* --chr : Path to a file that contains a list chromosome names of from the reference fasta, which should be included in the variant calling.
* --split : Length in BP of section in which the scaffolds of the reference are divided will be variant called separately.
* --poplist : Path to tab separated file with columns. Column 1 contains the the sample names and column two contains the name of the population the sample belongs to. Using this parameter initiates variant calling for each population and not all samples as one population. 
* -with-report : Nextflow parameter to create a html report about the computational resources the job needed. (Not to be confused with the html report that will be created from the variant calling statistics. That one will be in the outdir.)

**Splitting Scaffolds:**

The scaffolds are split into sections of the size given by the split parameter. This allows the parallelization of the variant calling step. According to the length chosen, each scaffold will be split into to so many sections and a variant call job will be started for each of them. Since the whole pipeline is per default restricted to start only 5 jobs at once and each variant calling job currently uses 100 GB of RAM, all of these settings need to be adapted to the analyzed data set. (RAM and queue size can be changed in the mfn.config file in the config folder.) Otherwise, the pipeline can easily get killed because it exceeds the memory limit, but if thats happening you can resume from were it stopped adding the "-resume" parameter.  

**Output:**

The output consists of one final vcf file holding all chromosomes and individuals. In addition to this there will be a Report.html file, which contains a summary from the variant statistics. The statistics files from which the report was created, are also part of the output and will be located in a folder called stats.

## Variant filtering

In addition to the variant calling pipeline, there is a small pipeline in the filter_variant folder for filtering the variants. It includes these steps:

1. Filtering of the variants using vcftools.
2. Calculating statistics about the file using vcftools and vcflib.
3. Creating a html report from the stats.

### Usage

```bash
nextflow run filter.nf --vcf /PATH/TO/VCF.vcf --outdir /PATH/TO/OUTPUT/DIRECTORY/ --frac 0.01

```

**Parameters:**

*  --vcf : Path to to the vcf file.
*  --outdir : Path to a directory were the output should be written into.
*  --frac : A fration to which the VCF file should be subsampled for the more precise statistics. Saves time for large VCF's.

In addition to those parameters, filter parameters can be set. The parameter defaults can be found in the nextflow.config file inside the filter_variants folder. I find it more convient to set them there, but they can also be set in the command. The parameters that can be set are the following.

*  --min_qual : Minimum variant quality filter.
*  --min_depth : Minimum read depth at variant position.
*  --max_depth : Maximum read depth at variant position.
*  --missingness : Fraction of not missing individuals per variant.
*  --min_maf : Minimum frequency of the minor allele.

Per default, indels are also removed in the filtering step.

**Output:**

The main output is the filtered vcf file and there will be the same stats and report files for the filtered file as were created in the variant calling pipeline. So maybe don't use the same ouput folder for the filter and variant calling pipelines, otherwise stats files might be overwritten.
