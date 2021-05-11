#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { nanopore } from './modules/nanopore.nf'

fast5 = Channel.fromPath(params.sequencing.fast5,type: 'dir')
fastq = params.sequencing.fastq == null ? '' : Channel.fromPath(params.sequencing.fastq, type: 'dir')
data = Channel.fromPath("${params.data}")

workflow {
  nanopore(fast5, fastq,data)
}
