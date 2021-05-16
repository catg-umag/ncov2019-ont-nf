nextflow.enable.dsl = 2


workflow Assembly {
  take:
    sample_names
    fast5_dir
    fastq_dirs
    sequencing_summary

   main:
    if (fastq_dirs == null) {
      basecalling(fast5_dir)

      demultiplexing(basecalling.out.all)
        | flatten
        | set { fastq_dirs }

      sequencing_summary = basecalling.out.sequencing_summary
    }

    fastq_dirs
      | map { [it.baseName, it] }
      | join(sample_names)
      | filtering

    articConsensus(
      fast5_dir,
      sequencing_summary,
      filtering.out.flat
    )

    articConsensus.out.consensus
      | map { it[1] }
      | collectFile(name: "all_consensus.fasta", newLine: true)
      | collect
      | set { all_consensus }

    downloadSnpEffDb()
    annotateVCFs(articConsensus.out.vcf, downloadSnpEffDb.out)

  emit:
    consensus = all_consensus
    bam = articConsensus.out.bam
    vcf = annotateVCFs.out
}


process basecalling {
  label 'guppy'
  label 'gpu'
  cpus 16
  
  input:
  path(fast5_dir)

  output:
  path('basecalled/'), emit: all
  path('basecalled/sequencing_summary.txt'), emit: sequencing_summary

  script:
  """
  guppy_basecaller \
    --input_path ${fast5_dir} \
    --save_path basecalled \
    --config ${params.guppy_basecalling_config} \
    --recursive \
    --device ${params.guppy_device} \
    --num_callers ${task.cpus} \
    --gpu_runners_per_device 4 \
    --chunks_per_runner 1024 \
    --chunk_size 2048
  """
}


process demultiplexing {
  label 'guppy'
  cpus 16

  input:
  path(fastq_dir)

  output:
  path('demultiplexed/barcode*')
  
  script:
  both_ends = params.guppy_demux_both_ends ? '--require_barcodes_both_ends' : ''
  """
  guppy_barcoder \
    --input_path ${fastq_dir}/pass \
    --save_path demultiplexed/ \
    --recursive \
    --barcode_kits "${params.guppy_barcodes}" \
    $both_ends \
    --worker_threads ${task.cpus}
  """
}


process filtering {
  tag "$sample"
  label 'artic'
  publishDir "${params.output_directory}/raw_data/", mode: 'copy', pattern: '*.fastq.gz'
  cpus 1

  input:
  tuple val(barcode), path(fastq_dir), val(sample)

  output:
  tuple val(sample), path("${sample}.fastq"), emit: flat
  tuple val(sample), path("${sample}.fastq.gz"), emit: compressed

  script:
  """
  artic guppyplex --min-length 400 --max-length 700 \
    --directory $fastq_dir --output ${sample}.fastq

  gzip -c ${sample}.fastq > ${sample}.fastq.gz
  """
}


process articConsensus {
  tag "$sample"
  label 'artic'
  publishDir "${params.output_directory}/artic/${sample}", mode: 'copy'
  cpus 4

  input:
  path(fast5_dir)
  path(sequencing_summary)
  tuple val(sample), path(fastq_file)

  output:
  tuple val(sample), path('*')
  tuple val(sample), path('*.consensus.fasta'), emit: consensus
  tuple val(sample), path('*.primertrimmed.rg.sorted.bam'), emit: bam
  tuple val(sample), path('*.pass.vcf.gz'), emit: vcf

  script:
  """
  artic minion --threads ${task.cpus} \
    --normalise ${params.artic_normalise} \
    --scheme-directory /opt/artic-ncov2019/primer_schemes/ \
    --fast5-directory $fast5_dir \
    --sequencing-summary $sequencing_summary \
    --read-file $fastq_file \
    ${params.artic_primer_scheme} $sample
  """
}


process downloadSnpEffDb {
  label 'snpeff'

  output:
  path('snpeff_data')

  script:
  """
  snpEff download -dataDir \$PWD/snpeff_data NC_045512.2 
  """
}


process annotateVCFs {
  tag "$sample"
  label 'snpeff'
  publishDir "${params.output_directory}/vcf", mode: 'copy'

  input:
  tuple val(sample), path(vcf)
  path('snpeff_data')

  output:
  tuple val(sample), path("${sample}.anno.vcf")

  script:
  """
  zcat $vcf \
  | sed 's/MN908947.3/NC_045512.2/g' \
  | snpEff -canon -noLof -no-downstream -no-upstream -noStats -noLog \
    -nodownload -dataDir \$PWD/snpeff_data NC_045512.2 \
    > ${sample}.anno.vcf
  """
}