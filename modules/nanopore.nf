nextflow.enable.dsl = 2

workflow nanopore {
   take:
  fast5
  fastq

   main:
    basecalling(fast5)
    reads_to_filter = params.sequencing.basecalling == true ? basecalling.out.pass.flatten() : fastq
    sample_names = Channel.from(params.samples)
    reads_to_filter = reads_to_filter.map { it -> [it.baseName, it] }.join(sample_names)
    changeSampleName(reads_to_filter)
    filtering(changeSampleName.out)
    pipeline(filtering.out)
    sorted_bam = pipeline.out.bam.flatten().filter { !(it.Name.contains('trimmed')) }
    sorted_bam \
    | (depth & coverage)
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
  params.sequencing.basecalling

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

process changeSampleName {
  cpus 2
  publishDir "${params.outdir}/basecalling", mode: 'copy'
  input:
    tuple val(bc_id), path(fastq), val(sample_name)
  output:
    tuple val(bc_id), path(sample_name)

  script:
  """
  mv ${fastq} ${sample_name}
  """
}
process filtering {
  label 'artic'
  publishDir "${params.outdir}/artic/filter", mode: 'copy'
  cpus 2
  input:
    tuple val(bc_id), path(fastq)
  output:
    tuple val(bc_id), path('filter*')
  script:
  """
  artic guppyplex --min-length 400 --max-length 700 --directory ${fastq} --prefix filter
  """
}

process pipeline {
  label 'artic'
  cpus 4
  publishDir "${params.outdir}/artic/pipeline/${fastq.baseName}", mode: 'copy'
  publishDir "${params.outdir}/artic/pipeline/consensus/", mode: 'copy', pattern: '*.consensus.fasta'
  publishDir "${params.outdir}/artic/pipeline/vcf/", mode: 'copy', pattern: '*.pass.vcf.gz'
  publishDir "${params.outdir}/artic/pipeline/bam/", mode: 'copy', pattern: '*.sorted.bam'

  input:
  tuple val(bc_id), path(fastq)
  output:
  path('*')
  path("*.sorted.bam") , emit: bam

  script:
  id = fastq.baseName - ~/filter_/
  fast5 = params.sequencing.basecalling == true ? "${params.sequencing.fast5}" : "${params.sequencing.fast5}/${bc_id}"
  """
artic minion --normalise 200 --threads 4 --scheme-directory ${params.sequencing.primers} --read-file ${fastq} --fast5-directory ${fast5} --sequencing-summary ${params.sequencing.summary} ${params.sequencing.primers_scheme} ${fastq.baseName}
  """
}

process depth {
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

process coverage {
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
