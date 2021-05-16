nextflow.enable.dsl = 2


workflow LineageAssesment {
  take:
    consensus
    vcf
    gisaid_clades

  main:
    consensus
      | (assignPangolin & assignNextstrain)
    
    vcf
      | map { it[1] }
      | collect
      | set { all_vcfs }
    assignGisaid(all_vcfs, gisaid_clades)

    mergeLineages(
      assignPangolin.out,
      assignGisaid.out,
      assignNextstrain.out
    )

  emit:
    mergeLineages.out
}


process assignPangolin {
  label 'pangolin'
  publishDir "${params.output_directory}/lineages", mode: 'copy'

  input:
  path(consensus)

  output:
  path('pangolin_lineages.csv')
  
  script:
  """
  pangolin $consensus --outfile pangolin_lineages.csv
  """
}


process assignGisaid {
  label 'pandas'
  publishDir "${params.output_directory}/lineages", mode: 'copy'

  input:
  path(vcfs)
  path('gisaid_clade_mutations.csv')

  output:
  path('gisaid_clades.csv')

  script:
  """
  assign_clades.py \
    --clade-mutation-list gisaid_clade_mutations.csv \
    --output gisaid_clades.csv \
    --input-vcfs $vcfs
  """
}


process assignNextstrain {
  label 'nextclade'
  publishDir "${params.output_directory}/lineages", mode: 'copy'

  input:
  path(consensus)

  output:
  path("nextstrain_lineages.csv")

  script:
  data_dir = '/opt/nextclade_data/sars-cov-2'
  """
  nextclade --input-fasta=$consensus \
    --input-root-seq=${data_dir}/reference.fasta \
    --genes=E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
    --input-gene-map=${data_dir}/genemap.gff \
    --input-tree=${data_dir}/tree.json \
    --input-qc-config=${data_dir}/qc.json \
    --output-csv=nextstrain_lineages.csv \
    --output-dir=.
  """
}


process mergeLineages {
  label 'pandas'

  input:
  path(pangolin_lineages)
  path(gisaid_lineages)
  path(nextstrain_lineages)

  output:
  path('all_lineages.csv')

  script:
  """
  merge_lineage_info.py \
    --pangolin $pangolin_lineages \
    --gisaid $gisaid_lineages \
    --nextstrain $nextstrain_lineages \
    --output all_lineages.csv
  """
}