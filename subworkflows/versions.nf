nextflow.enable.dsl=2


workflow GetSoftwareVersions {
  (getArticVersion
    & getSnpEffVersion
    & getSamtoolsVersion
    & getPangolinVersion
    & getNextcladeVersion
  )
    | mix
    | set { versions }

  versions
    | collectFile(
        storeDir: "${params.output_directory}",
        newLine: true,
        sort: true,
        seed: 'Software\tVersion'
      ) {
        ['software_versions.tsv', "${it[0]}\t${it[1]}"]
      }
}


process getArticVersion {
  label 'artic'

  output:
  tuple val('ARTIC'), stdout

  script:
  """
  artic --version | grep -Eo '([0-9]+\\.)+[0-9]+' | tr -d '\n'
  """
}


process getSnpEffVersion {
  label 'snpeff'

  output:
  tuple val('SnpEff'), stdout

  script:
  """
  snpEff -version | grep -Eo '([0-9]+\\.)+[0-9]+[a-z]*' | tr -d '\n'
  """
}


process getSamtoolsVersion {
  label 'samtools'

  output:
  tuple val('Samtools'), stdout

  script:
  """
  samtools --version | head -1 | cut -d ' ' -f2 | tr -d '\n'
  """
}

process getPangolinVersion {
  label 'pangolin'

  output:
  tuple val('Pangolin'), stdout

  script:
  """
  pangolin --version | grep -Eo '([0-9]+\\.)+[0-9]+' | tr -d '\n'
  """
}


process getNextcladeVersion {
  label 'nextclade'

  output:
  tuple val('Nextclade'), stdout

  script:
  """
  nextclade --version | tr -d '\n'
  """
}