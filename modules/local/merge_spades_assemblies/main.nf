process MERGE_SPADES_ASSEMBLIES {
    tag "${meta.id}"
    label 'process_single'

    container 'python:3.12.11-bookworm'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), val(assembly_metadata), path(contigs, arity: '1..*', stageAs: 'assemblies/*')
    path merge_script

    output:
    tuple val(meta), path('spades_coassembly.contigs.fa'), emit: contigs
    tuple val(meta), path('spades_coassembly.contig_provenance.tsv'), emit: provenance
    tuple val(meta), path('spades_coassembly.assemblies.json'), emit: manifest
    tuple val("${task.process}"), val('merge_spades_assemblies'), val('1.0.0'), emit: versions

    script:
    def files = contigs instanceof List ? contigs : [contigs]
    def manifest = assembly_metadata.withIndex().collect { item, index ->
        [group: item.group, group_id: item.group_id, coassembly_id: item.id,
         sample_ids: item.sample_ids, path: files[index].toString()]
    }
    def manifest_shell = groovy.json.JsonOutput.toJson(manifest).replace("'", "'\"'\"'")
    """
    printf '%s\\n' '${manifest_shell}' > spades_coassembly.assemblies.json
    python3 "${merge_script}" merge \\
        --manifest spades_coassembly.assemblies.json \\
        --output spades_coassembly.contigs.fa \\
        --provenance spades_coassembly.contig_provenance.tsv
    """
}
