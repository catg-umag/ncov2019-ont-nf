nextflow.enable.dsl = 2


workflow GetStatistics {
  take:
    bamfiles  // channel [name, BAM file]

  main:
    alignmentStats(bamfiles)

    bamfiles
      | fastQC
      | collect
      | multiQC

    alignmentStats.out.coverage
      | collectFile(keepHeader: true, skip: 1) {
        ["coverages.tsv", it[1].text.replaceAll('MN908947.3', it[0])] }
      | collect
      | set { all_stats }

  emit:
    coverage = all_stats
}


process alignmentStats {
  label 'samtools'
  tag { sample }
  publishDir "${params.output_directory}/qc/alignment_stats", mode: 'copy'

  input:
  tuple val(sample), path(bamfile)

  output:
  tuple val(sample), path("${sample}.depth"), emit: depth
  tuple val(sample), path("${sample}.coverage"), emit: coverage

  script:
  """
  samtools depth -a ${bamfile} -o ${sample}.depth
  samtools coverage ${bamfile} -o ${sample}.coverage
  """
}


process fastQC {
  label 'fastqc'
  tag { sample }
  publishDir "${params.output_directory}/qc/fastqc", mode: 'copy'
  cpus 1

  input:
  tuple val(sample), path(bamfile)

  output:
  path("fastqc_${sample}")

  script:
  """
  mkdir fastqc_${sample}
  fastqc ${bamfile} -o fastqc_${sample} -t ${task.cpus}
  """
}


process multiQC {
  label 'multiqc'
  publishDir "${params.output_directory}/qc/multiqc", mode: 'copy'
  
  input:
  path(fastqc_reports)

  output:
  tuple path('multiqc_data'), path('multiqc_report.html')

  script:
  """
  cp ${workflow.projectDir}/conf/multiqc_config.yaml .
  multiqc .
  """
}