#!/usr/bin/env nextflow
log.info """\
         Fancy nextflow variant caller
         ===================================
         Bam directory           : ${params.bamdir}
         Reference file          : ${params.ref}
         Output directory        : ${params.outdir}
         Chromosome list         : ${params.chr}
         BP chromosome split     : ${params.split}
         Stats subsample fraction: ${params.frac}
         Populations list        : ${params.poplist}
         High cov number of reads: ${params.highcov}
         High cov divide by      : ${params.div}
         """
         .stripIndent()

// holds the individual bam files
Channel
  .fromPath("${params.bamdir}*/*.bam")
  .set{bam_ch}

Channel
  .fromPath("${params.bamdir}*/*.bam")
  .set{samp_ch}


chr_l = file(params.chr).readLines()
chr_list = Channel.from(chr_l)
//Channel
  //  .fromList(chr_l)
  //  .set(chr_list)
//chr_list.view {print "$it"}

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
  tuple file(bam), file('*.bai') into ibam2_ch

  script:
  """
  samtools index $bam
  """
}


process CreateChrlist {

  input:
  file(allf) from ibam2_ch.collect()
  each chr from chr_list

  output:
  stdout into chrsplit_ch
  stdout ch

  script:
  def bams = allf.findAll{it =~ /bam_RG$/}
  """
  split=\$(python ${projectDir}/bin/chr_splitter.py -c $chr -r $params.ref -p $params.split)

  if [[ -z "\$split" ]]
  then
      echo "${chr}"
  else
      for item in \$split; do
        mreads=\$(for f in ${bams.join(' ')}; do samtools view -c \$f \$item ; done | awk '{ sum += \$1 } END { if (NR > 0) print sum / NR }')
        if [ \$( echo "\$mreads < ${params.highcov}" | bc ) -ne 0 ] ; then
         echo \$item
        else
         start=\$(echo \$item | sed 's/.*://' | sed 's/-.*//')
         end=\$(echo \$item |sed 's/.*-//')
         inc=\$(echo "${params.split} / ${params.div}" | bc)
         for i in \$(seq \$start \$inc \$end); do
           e=\$(echo "(\$i + \$inc) - 1" | bc)
           if ((\$e > \$end)); then
            if ((\$end != \$i)); then
              echo "${chr}:\$i-\$end"
            fi
           else
            echo "${chr}:\$i-\$e"
          fi
         done
        fi
       done
  fi

   """
}

//ch.view { print "$it" }

chrsplit_lines = chrsplit_ch.splitText()

if (params.poplist == null) {
  process CallVariants {

    label 'RAM_high'

    input:
    file(allf) from ibam_ch.collect()
    each chr from chrsplit_lines

    output:
    file('*') into chrvcf_ch

    script:
    def bams = allf.findAll{it =~ /bam_RG$/}
    //def bams = allf.collect { assert it ==~ pattern }
    """
    freebayes -f ${params.ref} -r ${chr.trim()} ${bams.join(' ')} > ${chr.trim()}.vcf
    """

  }
} else {
  process CallVariantsPop {

    label 'RAM_high'

    input:
    file(allf) from ibam_ch.collect()
    each chr from chrsplit_lines

    output:
    file('*') into chrvcf_ch

    script:
    def bams = allf.findAll{it =~ /bam_RG$/}
    //def bams = allf.collect { assert it ==~ pattern }
    """
    freebayes -f ${params.ref} -r ${chr.trim()} ${bams.join(' ')} --populations ${params.poplist} > ${chr.trim()}.vcf
    """
  }
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

process SortSamples {

  input:
  file(samp_list) from samp_ch.collect()
  each file(chr) from chrvcfgz_ch

  output:
  file('*.gz') into chrvcfgzso_ch

  script:
  samps = samp_list.baseName
  """
  tabix ${chr}
  bcftools view ${chr} -s ${samps.join(',')} -Oz > S_${chr}
  """
}

process IndexVCF {

  input:
  file(chr) from chrvcfgzso_ch

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
  file('*final.vcf.gz') into vcf_ch

  script:
  def vcfs = all_chr.findAll{it =~ /gz$/}
  """
  bcftools concat -a ${vcfs.join(' ')} > nfvacal_out.vcf
  grep "^#" nfvacal_out.vcf > nfvacal_final.vcf
  grep -v "^#" nfvacal_out.vcf| sort -T ./ -k1,1V -k2,2g >> nfvacal_final.vcf
  bgzip nfvacal_final.vcf
  """
}

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
