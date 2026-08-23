process PHYLOPHLAN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.container ?: 'quay.io/biocontainers/phylophlan:3.1.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(mags, arity: '1..*')
    path phylophlan_db
    path phylophlan_config

    output:
    tuple val(meta), path('*.phylophlan.tree.nwk'), emit: tree
    tuple val(meta), path('*.phylophlan.alignment.fasta'), emit: alignment
    tuple val(meta), path('*.phylophlan'), emit: results
    tuple val(meta), path('*.phylophlan.log'), emit: log
    tuple val("${task.process}"), val('phylophlan'), val('3.1.1'), emit: versions
    tuple val("${task.process}"), val('iqtree'), val('3.0.1'), emit: versions_iqtree

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: meta.id
    def mag_files = mags instanceof List ? mags : [mags]

    if (!mag_files) {
        error 'PHYLOPHLAN requires at least one MAG FASTA file'
    }

    def entries = mag_files.collect { mag ->
        def raw_id = mag.name.replaceFirst(/(?i)\.(fa|fna|fasta)(\.gz)?$/, '')
        def mag_id = raw_id.replaceAll(/[^A-Za-z0-9_.-]/, '_')
        if (!mag_id) {
            error "PHYLOPHLAN could not derive an identifier from '${mag.name}'"
        }
        [file: mag, id: mag_id]
    }

    if (entries.collect { entry -> entry.id }.toSet().size() != entries.size()) {
        error 'PHYLOPHLAN MAG identifiers are not unique after sanitization'
    }

    def stage_commands = entries.collect { entry ->
        if (entry.file.name.toLowerCase().endsWith('.gz')) {
            return "gzip -cd \"${entry.file}\" > \"genomes/${entry.id}.fna\""
        }
        "ln -s \"\$PWD/${entry.file}\" \"genomes/${entry.id}.fna\""
    }.join('\n')

    def database_name = phylophlan_db.name

    """
    set -euo pipefail

    mkdir -p genomes databases
    ${stage_commands}

    if [[ ! -d "${phylophlan_db}" ]]; then
        echo "PhyloPhlAn database input must be a directory" >&2
        exit 1
    fi

    ln -s "\$PWD/${phylophlan_db}" "databases/${database_name}"

    iqtree_version=\$(iqtree --version 2>&1)
    if ! grep -Eq '(^|[^0-9])3\.0\.1([^0-9]|\$)' <<< "\${iqtree_version}"; then
        echo "PhyloPhlAn requires IQ-TREE 3.0.1; configure a combined container or use Conda" >&2
        exit 1
    fi

    echo "Running PhyloPhlAn phylogenomics"

    phylophlan \
        -i genomes \
        -d "${database_name}" \
        -t a \
        -f "${phylophlan_config}" \
        --databases_folder databases \
        --genome_extension .fna \
        --diversity high \
        --accurate \
        --nproc ${task.cpus} \
        --output_folder "${prefix}.phylophlan" \
        --verbose \
        ${args} \
        2>&1 | tee "${prefix}.phylophlan.log"

    mapfile -t tree_files < <(find "${prefix}.phylophlan" -maxdepth 1 -type f -name '*.treefile')
    if [[ "\${#tree_files[@]}" -ne 1 ]]; then
        echo "Expected one IQ-TREE treefile, found \${#tree_files[@]}" >&2
        exit 1
    fi
    cp "\${tree_files[0]}" "${prefix}.phylophlan.tree.nwk"

    mapfile -t alignments < <(find "${prefix}.phylophlan" -maxdepth 1 -type f -name '*_concatenated.aln')
    if [[ "\${#alignments[@]}" -ne 1 ]]; then
        echo "Expected one PhyloPhlAn concatenated alignment, found \${#alignments[@]}" >&2
        exit 1
    fi
    cp "\${alignments[0]}" "${prefix}.phylophlan.alignment.fasta"
    """

    stub:
    def prefix    = task.ext.prefix ?: meta.id
    def mag_files = mags instanceof List ? mags : [mags]
    def entries   = mag_files.collect { mag ->
        def raw_id = mag.name.replaceFirst(/(?i)\.(fa|fna|fasta)(\.gz)?$/, '')
        [id: raw_id.replaceAll(/[^A-Za-z0-9_.-]/, '_')]
    }
    def newick = '(' + entries.collect { entry -> "${entry.id}:0.1" }.join(',') + ');'
    def alignment = entries.collect { entry ->
        ">${entry.id}\nMSTUBSEQUENCEMSTUBSEQUENCE"
    }.join('\n')

    """
    set -euo pipefail

    mkdir -p "${prefix}.phylophlan"

    printf '%s\n' '${alignment}' > "${prefix}.phylophlan/genomes_concatenated.aln"

    printf '%s\n' '${newick}' > "${prefix}.phylophlan/genomes.tre.treefile"
    cp "${prefix}.phylophlan/genomes_concatenated.aln" "${prefix}.phylophlan.alignment.fasta"
    cp "${prefix}.phylophlan/genomes.tre.treefile" "${prefix}.phylophlan.tree.nwk"

    printf 'PhyloPhlAn stub phylogenomics completed\n' > "${prefix}.phylophlan/phylophlan.log"
    printf 'PhyloPhlAn stub phylogenomics completed\n' > "${prefix}.phylophlan.log"
    """
}
