nextflow.enable.dsl = 2


workflow GenerateSummaries {
  take:
    samples_data      // single CSV file containing sample information
    consensus         // single FASTA file containing all consensuses
    vcfs              // channel [name, VCF file]
    coverage_stats    // single CSV file containg coverage statistics
    lineages          // single CSV file containing lineages
    epicov_template   // single XLS file containing template for EpiCov upload

  main:
    vcfs
      | map { it[1] }
      | collect
      | gatherVcfInfo

    getReferenceCoveredStats(consensus)

    mergeSampleData(
      samples_data,
      lineages,
      coverage_stats,
      getReferenceCoveredStats.out.collect()
    )

    generateExcelSummary(mergeSampleData.out, gatherVcfInfo.out)

    prepareSubmission(
      epicov_template,
      samples_data,
      mergeSampleData.out,
      consensus
    )
}


process getReferenceCoveredStats {
  label 'python'

  input:
  path(consensuses)

  output:
  path('covered_stats.csv')
  
  script:
  """
  get_covered_stats.py -i $consensuses -o covered_stats.csv
  """
}


process gatherVcfInfo {
  label 'python'
  publishDir "${params.output_directory}/summary", mode: 'copy'

  input:
  path(vcfs)

  output:
  path("variants_list${params.run_suffix}.csv")

  script:
  """
  gather_vcf_info.py \
    --output variants_list${params.run_suffix}.csv \
    --input-vcfs $vcfs
  """
}


process mergeSampleData {
  label 'python'
  publishDir "${params.output_directory}/summary", mode: 'copy'

  input:
  path(base)
  path(lineages)
  path(coverage)
  path(ref_coverage)

  output:
  path("sample_summary${params.run_suffix}.csv")

  script:
  """
  merge_sample_data.py \
    --base-data $base \
    --lineages $lineages \
    --coverage-stats $coverage \
    --ref-coverage-stats $ref_coverage \
    --output sample_summary${params.run_suffix}.csv
  """
}


process generateExcelSummary {
  label 'python'
  publishDir "${params.output_directory}/summary", mode: 'copy'

  input:
  path(samples)
  path(variants)
  
  output:
  path("summary${params.run_suffix}.xlsx")
  
  script:
  """
  generate_excel_summary.py \
    --sample-summary $samples \
    --variants-list $variants \
    --required-ref-coverage ${params.publish_minimum_completion} \
    --output summary${params.run_suffix}.xlsx
  """
}


process prepareSubmission {
  label 'python'
  publishDir "${params.output_directory}/gisaid_submission", mode: 'copy'

  input:
  path(template)
  path(samples_data)
  path(samples_summary)
  path(sequences)

  output:
  tuple path("EpiCov_Upload${params.run_suffix}.xls"), \
        path("sequences${params.run_suffix}.fasta")

  script:
  """
  generate_submission.py \
    --submission-template $template \
    --base-data $samples_data \
    --sample-summary $samples_summary \
    --input-sequences $sequences \
    --required-ref-coverage ${params.publish_minimum_completion} \
    --output-excel EpiCov_Upload${params.run_suffix}.xlsx \
    --output-sequences sequences${params.run_suffix}.fasta

  # convert xlsx to xls and delete xlsx
  unoconv -o EpiCov_Upload${params.run_suffix}.xls EpiCov_Upload${params.run_suffix}.xlsx
  rm EpiCov_Upload${params.run_suffix}.xlsx
  """
}