workflow Assembly {
  take:
    samples             // channel [sample_id, fastq_file]
    fast5_dir           // single directory containing FAST5 files
    sequencing_summary  // single TXT file

   main:
    samples
      | filtering

    articConsensus(
      fast5_dir,
      sequencing_summary,
      filtering.out,
      getPrimerSchemes()
    )

    articConsensus.out.consensus
      | map { it[1] }
      | collectFile(
          name: "all_consensus${params.run_suffix}.fasta",
          storeDir: "${params.output_directory}/summary",
          newLine: true)
      | collect
      | set { all_consensus }

    downloadSnpEffDb()
    annotateVCFs(articConsensus.out.vcf, downloadSnpEffDb.out)

  emit:
    consensus = all_consensus
    bam = articConsensus.out.bam
    vcf = annotateVCFs.out
}


process filtering {
  tag { sample }
  label 'fastp'
  publishDir "${params.output_directory}/raw_data/", mode: 'copy', pattern: '*.fastq.gz'
  cpus 1

  input:
  tuple val(sample), path(fastq_file)

  output:
  tuple val(sample), path("filtered/${sample}.fastq.gz")

  script:
  """
  mkdir filtered
  fastp \
    -i ${fastq_file} \
    -o filtered/${sample}.fastq.gz \
    --length_required 400 \
    --length_limit 700 \
    --disable_adapter_trimming \
    --disable_quality_filtering
  """
}


process getPrimerSchemes {
  label 'git'

  output:
  path('primer-schemes')

  script:
  """
  git clone https://github.com/artic-network/primer-schemes/ 
  """
}


process articConsensus {
  tag { sample }
  label 'artic'
  publishDir "${params.output_directory}/artic/${sample}", mode: 'copy'
  cpus 4
  memory 6.GB

  input:
  path(fast5_dir)
  path(sequencing_summary)
  tuple val(sample), path(fastq_file)
  path(primer_schemes)

  output:
  tuple val(sample), path('*')
  tuple val(sample), path('*.consensus.fasta'), emit: consensus
  tuple val(sample), path('*.primertrimmed.rg.sorted.bam'), emit: bam
  tuple val(sample), path('*.pass.vcf.gz'), emit: vcf

  script:
  variant_tool_opts = params.artic_use_nanopolish
    ? "--fast5-directory ${fast5_dir} --sequencing-summary ${sequencing_summary}"
    : "--medaka --medaka-model ${params.artic_medaka_model}"
  """
  artic minion \
    --threads ${task.cpus} \
    --normalise ${params.artic_normalise} \
    --scheme-directory ${primer_schemes} \
    ${variant_tool_opts} \
    --read-file ${fastq_file} \
    ${params.artic_primer_scheme} ${sample}
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
  tag { sample }
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