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
Channel
  .fromPath(params.sample_data)
  .splitCsv(header: true)
  .map { row -> [row.sample, pathCheck(row.fastq_file)] }
  .set { samples }

sequencing_summary = params.sequencing_summary != null
  ? pathCheck(params.sequencing_summary)
  : null

fast5_dir = pathCheck(params.fast5_directory) != null
  ? pathCheck(params.sequencing_summary)
  : null
epicov_template = pathCheck(params.gisaid_template)


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
    samples_data,
    Assembly.out.consensus,
    Assembly.out.vcf,
    GetStatistics.out.coverage,
    LineageAssesment.out.lineages,
    epicov_template
  )

  GetSoftwareVersions()
}