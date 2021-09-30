#!/usr/bin/env nextflow
log.info """\
         Fancy nextflow variant caller
         ===================================
         Bam directory           : ${params.bamdir}
         Reference file          : ${params.ref}
         Output directory        : ${params.outdir}
         Chromosome list         : ${params.chr}
         """
         .stripIndent()

// holds the individual bam files
Channel
  .fromPath("${params.bamdir}*/*.bam")
  .set{bam_ch}

chr_list = file(params.chr).readLines()

process AddRG {

  input:
  file(bam) from bam_ch

  output:
  file('*') into rgbam_ch

  script:
  """
  bamaddrg -b $bam -s ${bam.baseName} > ${bam}_RG
  """
}

process IndexBams {

  input:
  file(bam) from rgbam_ch

  output:
  tuple file(bam), file('*.bai') into ibam_ch

  script:
  """
  samtools index $bam
  """
}

process CallVariants {

  input:
  file(allf) from ibam_ch.collect()
  each chr from chr_list

  output:
  file('*') into chrvcf_ch

  script:
  def bams = allf.findAll{it =~ /bam_RG$/}
  //def bams = allf.collect { assert it ==~ pattern }
  """
  freebayes -f ${params.ref} -r $chr ${bams.join(' ')}  > ${chr}.vcf
  """


}

process CompressVCF {

  input:
  file(chr) from chrvcf_ch

  output:
  file('*') into chrvcfgz_ch

  script:
  """
  bgzip -c $chr > ${chr}.gz
  """
}

process IndexVCF {

  input:
  file(chr) from chrvcfgz_ch

  output:
  tuple file(chr), file('*') into chrvcfgzi_ch

  script:
  """
  bcftools index $chr
  """
}

process MergeVCF {

  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(all_chr) from chrvcfgzi_ch.collect()

  output:
  file('*out.vcf.gz*') into vcf_ch

  script:
  def vcfs = all_chr.findAll{it =~ /gz$/}
  """
  bcftools concat ${vcfs.join(' ')} | bgzip -c > nfvacal_out.vcf.gz
  bcftools index nfvacal_out.vcf.gz
  """
}
