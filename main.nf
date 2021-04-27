#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { nanopore } from './modules/nanopore.nf'

fast5 = Channel.fromPath(params.sequencing.fast5,type: 'dir')
fastq= fast5
if(!params.sequencing.basecalling)
  fastq = Channel.fromPath(params.sequencing.fastq, type: 'dir')

workflow{
  nanopore(fast5,fastq)
}
