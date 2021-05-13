nextflow.enable.dsl=2


workflow GenerateSummaries {
  take:
    samples_data
    consensus
    coverage_stats
    lineages
    epicov_template

  main:
    getReferenceCoveredStats(consensus)

    mergeData(
      samples_data,
      lineages,
      coverage_stats,
      getReferenceCoveredStats.out.collect())
}


process getReferenceCoveredStats {
  label 'biopython'

  input:
  path(consensuses)

  output:
  path('covered_stats.csv')
  
  script:
  mincov = params.publish_minimum_coverage
  """
  get_covered_stats.py -i $consensuses -o covered_stats.csv
  """
}


process mergeData {
  label 'pandas'
  publishDir "${params.output_directory}", mode: 'copy'

  input:
  path(base)
  path(lineages)
  path(coverage)
  path(ref_coverage)

  output:
  path('summary.csv')

  script:
  """
  merge_sample_data.py \
    --base-data $base \
    --lineages $lineages \
    --coverage-stats $coverage \
    --ref-coverage-stats $ref_coverage \
    --output summary.csv
  """
}