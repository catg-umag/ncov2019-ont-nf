#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { addDefaultParamValues } from './lib/groovy/utils.gvy'

addDefaultParamValues(params, "${workflow.projectDir}/params.default.yml")
params.run_suffix = params.run_id != null ? "_${params.run_id}" : ''

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

if (params.fastq_directory != null) {
  Channel
    .fromPath("${params.fasq_directory}/barcode*", type: 'dir')
    .set { fastq_dirs }
} else {
  fastq_dirs = null
}

sequencing_summary = params.sequencing_summary != null
  ? file(params.sequencing_summary)
  : null

samples_data = file(params.sample_data)
fast5_dir = file(params.fast5_directory)
gisaid_clades = file(params.gisaid_clades)
epicov_template = file(params.gisaid_template)


workflow {
  Assembly(
    sample_names,
    fast5_dir,
    fastq_dirs,
    sequencing_summary
  )

  GetStatistics(Assembly.out.bam)

  LineageAssesment(
    Assembly.out.consensus,
    Assembly.out.vcf,
    gisaid_clades
  )

  GenerateSummaries(
    samples_data,
    Assembly.out.consensus,
    Assembly.out.vcf,
    GetStatistics.out.coverage,
    LineageAssesment.out,
    epicov_template
  )
}