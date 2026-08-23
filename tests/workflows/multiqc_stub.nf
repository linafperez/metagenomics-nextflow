#!/usr/bin/env nextflow

include { GLOBAL_PROCESSING_EVALUATION } from '../../subworkflows/local/global_processing_evaluation/main'

workflow {
    def pipeline_root = "${projectDir}/../.."

    ch_reports = channel.of(
        tuple(
            [id: 'arbitrary_reports'],
            [
                file("${pipeline_root}/assets/phylophlan_iqtree.cfg", checkIfExists: true),
                file("${pipeline_root}/assets/multiqc_config.yml", checkIfExists: true)
            ]
        )
    )

    ch_multiqc_config = channel.value(
        file("${pipeline_root}/assets/multiqc_config.yml", checkIfExists: true)
    )

    GLOBAL_PROCESSING_EVALUATION(ch_reports, ch_multiqc_config)
}
