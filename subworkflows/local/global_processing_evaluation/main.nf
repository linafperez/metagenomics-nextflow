#!/usr/bin/env nextflow

include { MULTIQC } from '../../../modules/core/multiqc/main'

def collect_report_paths(item) {
    if (item == null || item instanceof Map) {
        return []
    }
    if (item instanceof java.nio.file.Path) {
        return [item]
    }
    if (item instanceof File) {
        return [item.toPath()]
    }
    if (item instanceof Collection) {
        return item.collectMany { nested -> collect_report_paths(nested) }
    }
    return []
}

workflow GLOBAL_PROCESSING_EVALUATION {
    take:
    ch_reports
    ch_multiqc_config

    main:
    ch_report_files = ch_reports
        .flatMap { item -> collect_report_paths(item) }
        .unique()

    ch_multiqc_input = ch_report_files
        .collect()
        .map { reports ->
            if (!reports) {
                error 'GLOBAL_PROCESSING_EVALUATION requires at least one report file or directory'
            }
            tuple(
                [id: 'global_processing_evaluation', stage: 'global_processing_evaluation'],
                reports
            )
        }

    MULTIQC(ch_multiqc_input, ch_multiqc_config)

    emit:
    report   = MULTIQC.out.report
    data     = MULTIQC.out.data
    logs     = MULTIQC.out.log
    versions = MULTIQC.out.versions
}
