process INTEGRATE_ANNOTATIONS {
    tag "${meta.mag_id ?: meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.container ?: 'python:3.12.11-slim-bookworm'}"

    input:
    tuple val(meta), path(proteins), path(gff), path(eggnog_annotations), path(interproscan_tsv)
    path annotation_script

    output:
    tuple val(meta), path('*.functional_annotations.tsv'), emit: annotations
    tuple val(meta), path('*.functional_annotations.summary.json'), emit: summary
    tuple val(meta), path('*.annotation_integration.log'), emit: log
    tuple val("${task.process}"), val('integrate_functional_annotations'), val('1.0.0'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: (meta.mag_id ?: meta.id)
    def mag_id = meta.mag_id ?: meta.id

    """
    set -euo pipefail

    python3 "${annotation_script}" \
        --proteins "${proteins}" \
        --gff "${gff}" \
        --eggnog "${eggnog_annotations}" \
        --interpro "${interproscan_tsv}" \
        --mag-id "${mag_id}" \
        --output "${prefix}.functional_annotations.tsv" \
        --summary "${prefix}.functional_annotations.summary.json" \
        ${args} \
        2>&1 | tee "${prefix}.annotation_integration.log"
    """

    stub:
    def prefix = task.ext.prefix ?: (meta.mag_id ?: meta.id)
    def mag_id = meta.mag_id ?: meta.id

    """
    set -euo pipefail

    protein_id=\$(awk '/^>/{sub(/^>/, ""); split(\$0, fields, /[[:space:]]+/); print fields[1]; exit}' "${proteins}")
    protein_id=\${protein_id:-stub_protein_1}

    printf '%b\n' 'mag_id\tprotein_id\tgene_id\tsequence_length_aa\tcontig\tfeature_start\tfeature_end\tstrand\tpreferred_name\tdescription\tseed_ortholog\teggnog_ogs\tmax_annotation_level\tcog_category\tgo_terms\tec_numbers\tkegg_ko\tkegg_pathways\tkegg_modules\tkegg_reactions\tcazy_families\tpfam_accessions\tinterpro_accessions\tinterpro_descriptions\tinterpro_member_databases\tinterpro_signature_accessions\tinterpro_signature_descriptions\tinterpro_pathways\tannotation_sources' \
        > "${prefix}.functional_annotations.tsv"
    printf '%s\t%s\t%s\t20\tstub_contig\t1\t63\t+\tstubA\tStub ribosomal protein\tstub.ortholog\tCOG0001@2|Bacteria\t2|Bacteria\tJ\tGO:0003674\t1.1.1.1\tko:K00001\tmap00010\tM00001\tR00001\tGH1\tPF00001\tIPR000001\tStub InterPro entry\tPfam\tPF00001\tStub protein domain\tKEGG:map00010\teggNOG-mapper;InterProScan\n' \
        '${mag_id}' "\${protein_id}" "\${protein_id}" >> "${prefix}.functional_annotations.tsv"

    printf '{"eggnog_annotated_proteins":1,"interpro_annotated_proteins":1,"mag_id":"%s","protein_count":1,"proteins_with_cazy":1,"proteins_with_ec":1,"proteins_with_go":1,"proteins_with_interpro":1,"proteins_with_kegg":1,"proteins_with_pfam":1}\n' \
        '${mag_id}' > "${prefix}.functional_annotations.summary.json"

    printf 'Functional annotation integration stub completed\n' > "${prefix}.annotation_integration.log"
    """
}
