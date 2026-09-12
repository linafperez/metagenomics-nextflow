#!/usr/bin/env nextflow

include { METAGENOMICS } from './workflows/metagenomics'

workflow {
    METAGENOMICS()
}
