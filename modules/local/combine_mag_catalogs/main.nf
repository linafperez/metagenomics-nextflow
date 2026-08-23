process COMBINE_MAG_CATALOGS {
    tag "${meta.id}"
    label 'process_single'

    container 'python:3.12.11-slim-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta),
        path(megahit_mags, arity: '1..*', stageAs: 'megahit/*'),
        path(megahit_quality, stageAs: 'megahit_quality/*'),
        path(spades_mags, arity: '1..*', stageAs: 'spades/*'),
        path(spades_quality, stageAs: 'spades_quality/*')
    path combine_script

    output:
    tuple val(meta), path('combined_catalog/*.fa', arity: '2..*'), emit: mags
    tuple val(meta), path('combined_catalog.provenance.tsv'), emit: provenance
    tuple val(meta), path('combined_catalog.quality.tsv'), emit: quality
    tuple val(meta), path('combined_catalog.genomeInfo.csv'), emit: genome_info
    tuple val("${task.process}"), val('combine_mag_catalogs'), val('1.0.0'), emit: versions

    script:
    def megahitFiles = megahit_mags instanceof List ? megahit_mags : [megahit_mags]
    def spadesFiles = spades_mags instanceof List ? spades_mags : [spades_mags]
    def megahitArgs = megahitFiles.collect { mag -> "\"${mag}\"" }.join(' ')
    def spadesArgs = spadesFiles.collect { mag -> "\"${mag}\"" }.join(' ')
    def command = """
    python3 "${combine_script}" \\
        --megahit-bins ${megahitArgs} \\
        --spades-bins ${spadesArgs} \\
        --megahit-quality "${megahit_quality}" \\
        --spades-quality "${spades_quality}" \\
        --output-dir combined_catalog \\
        --provenance combined_catalog.provenance.tsv \\
        --quality-table combined_catalog.quality.tsv \\
        --genome-info combined_catalog.genomeInfo.csv
    """
    command

    stub:
    def megahitFiles = megahit_mags instanceof List ? megahit_mags : [megahit_mags]
    def spadesFiles = spades_mags instanceof List ? spades_mags : [spades_mags]
    def megahitArgs = megahitFiles.collect { mag -> "\"${mag}\"" }.join(' ')
    def spadesArgs = spadesFiles.collect { mag -> "\"${mag}\"" }.join(' ')
    """
    python3 "${combine_script}" \\
        --megahit-bins ${megahitArgs} \\
        --spades-bins ${spadesArgs} \\
        --megahit-quality "${megahit_quality}" \\
        --spades-quality "${spades_quality}" \\
        --output-dir combined_catalog \\
        --provenance combined_catalog.provenance.tsv \\
        --quality-table combined_catalog.quality.tsv \\
        --genome-info combined_catalog.genomeInfo.csv
    """
}
