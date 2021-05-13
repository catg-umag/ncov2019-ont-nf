#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { addDefaultParamValues } from './lib/groovy/utils.gvy'

addDefaultParamValues(params, "${workflow.projectDir}/params.default.yml")

// load subworkflows after modifying parameters
include { Assembly } from './modules/assembly.nf'
include { GetStatistics } from './modules/statistics.nf'
include { LineageAssesment } from './modules/lineages.nf'
include { GenerateSummaries } from './modules/summaries.nf'


Channel
  .fromPath(params.sample_data)
  .splitCsv(header: true)
  .map { row -> [row.barcode, row.sample] }
  .set { sample_names }

fast5_dir = file(params.fast5_directory)
fastq_dirs = params.fastq_directory != null ? Channel.fromPath(params.fastq_directory, type: 'dir') : null
sequencing_summary = params.sequencing_summary != null ? file(params.sequencing_summary) : null
primer_schemes_dir = file(params.artic_primer_schemes_directory)
gisaid_clades = file(params.gisaid_clades)
samples_data = file(params.sample_data)
epicov_template = file(params.gisaid_template)


workflow {
  Assembly(
    sample_names,
    fast5_dir,
    fastq_dirs,
    sequencing_summary,
    primer_schemes_dir)

  GetStatistics(Assembly.out.bam)

  LineageAssesment(
    Assembly.out.consensus,
    Assembly.out.vcf,
    gisaid_clades)

  GenerateSummaries(
    samples_data,
    Assembly.out.consensus,
    GetStatistics.out.coverage,
    LineageAssesment.out,
    epicov_template)
}