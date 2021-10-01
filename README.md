# Nfvacal - nextflow variant calling pipeline using freebayes

 This pipeline calls variants from bam files of multiple individuals using freebayes. It is adapted for usage on the MfN cluster. 

Steps of the pipeline are: 

1. Adding RG individual tags to the bam file to later be able to identify different individuals in the final cvf file.
2. Index bam files./
3. Calling of variants for all individuals simultaneously, but in parallel for the different chromosomes of the reference genome using freebayes. 
4. Compressing the resulting vcf files. 
5. Merging of all chromosome vcf's into one final file. 
6. Calculating statistics about the file using vcftools and vcflib.
7. Creating a html report from the stats. 



## Usage

```bash
nextflow run vacal.nf --bamdir /PATH/TO/BAMS/ --ref /PATH/TO/REFERENCE/file.fasta --outdir /PATH/TO/DIR/TO/PUT/OUTPUT/ --chr /PATH/TO/CHROMOSOME/file.list -with-report /PATH/TO/REPORT/file.html
```

**Parameters:**

* --bamdir : Path to directory in with sub folders for each individual. Each sub folder must hold one bam file. 
* --ref : Path to the reference genome assembly fasta file.
* --outdir : Path to directory were the output of the pipeline will be written into. 
* --chr : Path to a file that contains a list chromosome names of from the reference fasta, which should be included in the variant calling. 
* -with-report : Nextflow parameter to create a html report about the computational resources the job needed. (Not to be confused with the html report that will be created from the variant calling statistics. That one will be in the outdir.)



**Output:**

The output consists of one final vcf file holding all chromosomes and individuals. In addition to this there will be a Report.html file, which contains a summary from the variant statistics. The statistics files from which the report was created, are also part of the output and will be located in a folder called stats. 