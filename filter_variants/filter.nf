#!/usr/bin/env nextflow
log.info """\
         Fancy nextflow variant filterer
         ===================================
         VCF file                : ${params.vcf}
         Out directory           : ${params.outdir}
         Stats subsample fraction: ${params.frac}
         ===================================
         Filter:
         -----------------------------------
         Minimum quality         : ${params.min_qual}
         Minimum depth           : ${params.min_depth}
         Maximum depth           : ${params.max_depth}
         Missingness             : ${params.missingness}
         Minimum MAF             : ${params.min_maf}
         """
         .stripIndent()

Channel
  .fromPath("${params.vcf}")
  .set{vcf_ch}


process Stats {
  publishDir "${params.outdir}/stats/", mode: 'copy'

  input:
  file(vcf) from vcf_ch

  output:
  tuple file('*.frq'), file('*.het'), file('*.idepth'), file('*.imiss'), file('*.ldepth.mean'), file('*.lmiss'), file('*.lqual'), file('*.txt') into stats_ch

  script:
  """
  bcftools index -f $vcf
  tabix $vcf
  vcfstats $vcf > vcfstats.txt
  vcfrandomsample -r ${params.frac}  $vcf > ${vcf.baseName}_subset.vcf
  bgzip ${vcf.baseName}_subset.vcf
  bcftools index -f ${vcf.baseName}_subset.vcf.gz
  tabix ${vcf.baseName}_subset.vcf.gz
  vcftools --gzvcf ${vcf.baseName}_subset.vcf.gz --freq2 --max-alleles 2
  vcftools --gzvcf ${vcf.baseName}_subset.vcf.gz --depth
  vcftools --gzvcf ${vcf.baseName}_subset.vcf.gz --site-mean-depth
  vcftools --gzvcf ${vcf.baseName}_subset.vcf.gz --site-quality
  vcftools --gzvcf ${vcf.baseName}_subset.vcf.gz --missing-indv
  vcftools --gzvcf ${vcf.baseName}_subset.vcf.gz --missing-site
  vcftools --gzvcf ${vcf.baseName}_subset.vcf.gz --het

  """
}

process Report {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple file(frq), file(het), file(idepth), file(imiss), file(ldepth), file(lmiss), file(lqual), file(vcfstats) from stats_ch

  output:
  file('*.html') into report_ch

  script:
  """
  Rscript -e "rmarkdown::render('${projectDir}/bin/Report.Rmd', params=list(vcfstats='\$PWD/${vcfstats}',lqual='\$PWD/${lqual}',idepth='\$PWD/${idepth}',ldepth='\$PWD/${ldepth}',imiss='\$PWD/${imiss}',lmiss='\$PWD/${lmiss}',het='\$PWD/${het}',frq='\$PWD/${frq}'), output_dir='\$PWD')"
  """

}
