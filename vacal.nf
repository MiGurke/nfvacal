#!/usr/bin/env nextflow
log.info """\
         Fancy nextflow variant caller
         ===================================
         Bam directory           : ${params.bamdir}
         Reference file          : ${params.ref}
         Output directory        : ${params.outdir}
         Chromosome list         : ${params.chr}
         BP chromosome split     : ${params.split}
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

process CreateChrlist {
  publishDir "${workDir}", mode: 'copy'

  input:
  val(chr) from chr_list

  output:
  stdout into chrsplit_ch

  script:
  """
#!/usr/bin/env python
from Bio import SeqIO

chr = "$chr"
ref = SeqIO.parse("$params.ref","fasta")
piece = $params.split

for rec in ref:
    id = rec.id
    if id == chr:
        nparts = len(rec.seq)/piece
        if nparts < 1:
            line = "$chr"
            print(line)
        else:
            for i in range(1,int(nparts) + 1):
                if i == 1:
                    start = 0
                    end = piece * i
                elif i == int(nparts):
                    start = piece * i + 1
                    end = len(rec.seq)
                else:
                    start = piece * (i-1) + 1
                    end = piece * i
                line = "$chr"+":"+str(start)+"-"+str(end)
                print(line)

  """
}


chrsplit_lines = chrsplit_ch.splitText()

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
  freebayes -f ${params.ref} -r ${chr.trim()} ${bams.join(' ')}  > ${chr.trim()}.vcf
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
  file('*final.vcf.gz') into vcf_ch

  script:
  def vcfs = all_chr.findAll{it =~ /gz$/}
  """
  bcftools concat -a ${vcfs.join(' ')} > nfvacal_out.vcf
  grep "^#" nfvacal_out.vcf > nfvacal_final.vcf
  grep -v "^#" nfvacal_out.vcf| sort -k1,1V -k2,2g >> nfvacal_final.vcf
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
  vcfrandomsample -r 0.01  $vcf > ${vcf.baseName}_subset.vcf
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
