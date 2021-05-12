nextflow.enable.dsl = 2

workflow nanopore {
   take:
    fast5
    fastq
    data
    gisaid_clades,
    gisaid_template

   main:
    sample_names = data.splitCsv().map { it -> [it[1], it[0]] }

    basecalling(fast5)
    reads_to_filter = params.sequencing.fastq == null ? basecalling.out.pass.flatten() : fastq
    reads_to_filter
  | map { it -> [it.baseName, it] }
  | join(sample_names)
  | filtering
  pipeline(filtering.out.filter)

  pipeline.out.bam
   | flatten
   | filter { !(it.Name.contains('trimmed')) }
   | (calculateDepth & coverageStats)

  assign_pangolin(pipeline.out.consensus.collect())
  pangolin = assign_pangolin.out.short_pangolin.splitCsv()

  assign_clade(pipeline.out.vcf)
  clades = assign_clade.out.map { name, clade -> [name, clade.readLines()[0]] }

  report_data = clades.join(pangolin).map { name, clade, pang -> [name.replaceAll(/filter_/, ''), clade, pang] }
  report_data = data.splitCsv().join(report_data)

  name = coverageStats.out.sample_names.map { it -> it }.collect()
  cov = coverageStats.out.cov.map {  it.readLines()[1].tokenize("\t")[6] }.collect()

  fillData(name, cov, gisaid_template)
}

process basecalling {
  label 'guppy'
  label 'gpu'
  cpus 16
  input:
    path(fast5)
  output:
  path("basecalled_${params.sequencing.id}")
  path("basecalled_${params.sequencing.id}/pass/barcode*"), emit: pass

  when:
  params.sequencing.fastq == null

  script:
  """
    guppy_basecaller \
    --input_path ${fast5} \
    --save_path basecalled_${params.sequencing.id} \
    --config ${params.sequencing.guppy_config} \
    --recursive \
    --device ${params.sequencing.guppy_device} \
    --num_callers ${task.cpus} \
    --gpu_runners_per_device 4 \
    --chunks_per_runner 1664 \
    --barcode_kits ${params.sequencing.barcodes}
  """
}

process filtering {
  tag "${sample_name}"
  label 'artic'
  publishDir "${params.outdir}/basecalling", mode: 'copy', pattern: "${sample_name}"
  publishDir "${params.outdir}/artic/filter", mode: 'copy'
  cpus 2
  input:
    tuple val(bc_id), path("${sample_name}"), val(sample_name)
  output:
    tuple val(sample_name), path('filter*'), val(bc_id), emit: filter
    path("${sample_name}")

  script:
  """
  artic guppyplex --min-length ${params.len_min_amplicon} --max-length ${params.len_max_amplicon} --directory ${sample_name} --prefix filter
  """
}

process pipeline {
  tag "${sample_name}"
  label 'artic'
  cpus 4
  publishDir "${params.outdir}/artic/pipeline/${fastq.baseName}", mode: 'copy'
  publishDir "${params.outdir}/artic/pipeline/consensus/", mode: 'copy', pattern: '*.consensus.fasta'
  publishDir "${params.outdir}/artic/pipeline/vcf/", mode: 'copy', pattern: '*.pass.vcf.gz'
  publishDir "${params.outdir}/artic/pipeline/bam/", mode: 'copy', pattern: '*.sorted.bam'

  input:
  tuple val(sample_name), path(fastq), val(bc_id)

  output:
  path('*')
  path("*.sorted.bam") , emit: bam
  path("*.consensus.fasta"), emit: consensus
  path("*.pass.vcf.gz"), emit: vcf

  script:
  """
    artic minion --normalise 200 --threads 4 --scheme-directory ${params.sequencing.primers} --read-file ${fastq} --fast5-directory ${params.sequencing.fast5} --sequencing-summary ${params.sequencing.summary} ${params.sequencing.primers_scheme} ${fastq.baseName}
  """
}

process calculateDepth {
  tag "${bam.baseName}"
  label 'samtools'
  cpus 2
  publishDir "${params.outdir}/depths",mode: 'copy'
  input:
  path(bam)

  output:
  path("${bam.simpleName}.depth")
  script:
  """
  samtools depth ${bam} -o ${bam.simpleName}.depth
  """
}

process coverageStats {
  tag "${bam.baseName}"
  label 'samtools'
  cpus 2
  publishDir "${params.outdir}/coverage",mode: 'copy'
  input:
  path(bam)

  output:
  path("${bam.simpleName}.coverage")
  script:
  """
  samtools coverage ${bam} -o ${bam.simpleName}.coverage
  """
}

process assign_pangolin {
  cpus 2
  label 'pangolin'
  publishDir "${params.outdir}/", mode: 'copy', pattern: 'all_consensus.fasta'
  publishDir "${params.outdir}/pangolin/", mode: 'copy', pattern: '*.csv'

  input:
  path(all_consensus)

  output:
  path('all_consensus.fasta')
  path('*.csv')
  path ("final_report_short.csv"), emit: short_pangolin
  script:
  """
  cat ${all_consensus} > all_consensus.fasta
  pangolin all_consensus.fasta --outfile lineages.csv
  cat lineages.csv | cut -d ',' -f1,2 > final_report_short.csv
  sed -i 's#/ARTIC/nanopolish_MN908947.3##g' final_report_short.csv
  """
}

process assign_clade {
  cpus 2

  label 'pandas'
  input:
  path(vcf)

  output:
  tuple val("${vcf.simpleName}"), path('clade.txt')
  script:
  """
  assign_clade.py ${vcf} ${params.gisaid_clades} > clade.txt
  """
}

// feo feo muy feo xD
process report {
  cpus 2
  label 'pandas'
  publishDir "${params.outdir}/", mode: 'copy'
  echo true
 input:
  tuple val(sample_name), val(barcode), val(name), val(clade), val(lineage) //colection_date,sex,age,location

 output:
  path('final_reportt.csv')

 script:
 """
 echo ${name},${clade},${lineage} >> ${params.report}
 cp ${params.report} final_reportt.csv
 """
}

process fillData {
  cpus 2
  label 'openpyxl'
  publishDir "${params.outdir}/excel_upload",mode: 'copy'
  echo true
  input:
  val(cov)
  val(name)
  path(gisaid_template)

  output:
  path('Submitted.xlsx')

  script:
  """
  fill_excel.py ${name} ${cov} ${params.gisaid_template}
  """
}
