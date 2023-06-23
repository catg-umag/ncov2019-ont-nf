#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { addDefaultParamValues; validateParameters; pathCheck } from './lib/groovy/utils.gvy'

// pre-process parameters (defaults and validations)
addDefaultParamValues(params, "${workflow.projectDir}/params.default.yml")
validateParameters(params, "${workflow.projectDir}/conf/parameter_config.yml")
params.run_suffix = params.run_id != null ? "_${params.run_id}" : ''

// load subworkflows after modifying parameters
include { Assembly } from './subworkflows/assembly.nf'
include { GetStatistics } from './subworkflows/statistics.nf'
include { LineageAssesment } from './subworkflows/lineages.nf'
include { GenerateSummaries } from './subworkflows/summaries.nf'
include { GetSoftwareVersions } from './subworkflows/versions.nf'

// transforms parameters in channels and variables
channel
  .fromPath(params.sample_data)
  .splitCsv(header: true)
  .map { row -> [row.sample, pathCheck(row.fastq_file)] }
  .set { samples }

sequencing_summary = params.sequencing_summary != null
  ? pathCheck(params.sequencing_summary)
  : file('NOFILE_SEQSUMM')

fast5_dir = params.fast5_directory != null
  ? pathCheck(params.sequencing_summary)
  : file('NOFILE_FAST5')
epicov_template = file(params.gisaid_template)

samples_file = file(params.sample_data)


workflow {
  Assembly(
    samples,
    fast5_dir,
    sequencing_summary
  )

  GetStatistics(Assembly.out.bam)

  LineageAssesment(
    Assembly.out.consensus,
  )

  GenerateSummaries(
    samples_file,
    Assembly.out.consensus,
    Assembly.out.vcf,
    GetStatistics.out.coverage,
    LineageAssesment.out.lineages,
    epicov_template
  )

  GetSoftwareVersions()
}