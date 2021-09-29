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
  tuple file(bam) into ibam_ch

  script:
  """
  samtools index $bam
  """
}

process CallVariants {
  echo true

  input:
  file(all_bams) from ibam_ch.collect()
  each chr from chr_list

  script:
  """
  freebayes -f ${params.ref} -r $chr $all_bams > ${params.outdir}/${chr}.vcf
  """


}
