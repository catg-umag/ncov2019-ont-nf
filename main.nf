#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { addDefaultParamValues; validateParameters } from './lib/groovy/utils.gvy'

// pre-process parameters (defaults and validations)
addDefaultParamValues(params, "${workflow.projectDir}/params.default.yml")
validateParameters(params, "${workflow.projectDir}/conf/parameter_config.yml")
params.run_suffix = params.run_id != null ? "_${params.run_id}" : ''

// load subworkflows after modifying parameters
include { Assembly } from './subworkflows/assembly.nf'
include { GetStatistics } from './subworkflows/statistics.nf'
include { LineageAssesment } from './subworkflows/lineages.nf'
include { GenerateSummaries } from './subworkflows/summaries.nf'

// transforms parameters in channels and variables
Channel
  .fromPath(params.sample_data)
  .splitCsv(header: true)
  .map { row -> [row.barcode, row.sample] }
  .set { sample_names }

fastq_dirs = params.fastq_directory != null
  ? Channel.fromPath("${params.fastq_directory}/barcode*", type: 'dir')
  : null

sequencing_summary = params.sequencing_summary != null
  ? file(params.sequencing_summary)
  : null

samples_data = file(params.sample_data)
fast5_dir = file(params.fast5_directory)
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