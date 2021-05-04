#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { nanopore } from './modules/nanopore.nf'

sample_names = Channel.from(params.samples)
fast5 = Channel.fromPath(params.sequencing.fast5,type: 'dir')
fastq = params.sequencing.basecalling ? '' : Channel.fromPath(params.sequencing.fastq, type: 'dir')

workflow {
  nanopore(fast5, fastq, sample_names)
}
