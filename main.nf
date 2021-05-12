#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { nanopore } from './modules/nanopore.nf'

fast5 = Channel.fromPath(params.sequencing.fast5,type: 'dir')
fastq = params.sequencing.fastq == null ? '' : Channel.fromPath(params.sequencing.fastq, type: 'dir')
clades_template = Channel.fromPath("${params.gisaid_clades}")
data = Channel.fromPath("${params.samples_data}")
gisaid_template = Channel.fromPath("${params.gisaid_template}")

workflow {
  nanopore(fast5, fastq, data, clades_template,gisaid_template)
}
