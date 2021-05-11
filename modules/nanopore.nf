nextflow.enable.dsl = 2

workflow nanopore {
   take:
    fast5
    fastq
    data

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
