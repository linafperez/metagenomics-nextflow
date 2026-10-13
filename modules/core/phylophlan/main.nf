process PHYLOPHLAN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.container ?: 'metagenomics/phylophlan-iqtree:3.1.1-3.0.1'}"

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

    NF_CPUS="${task.cpus}"
    NF_PREFIX="${prefix}"
    NF_DATABASE="${phylophlan_db}"
    NF_DATABASE_NAME="${database_name}"
    NF_CONFIG="${phylophlan_config}"

    mkdir -p genomes databases
    ${stage_commands}

    if [[ ! -d "\${NF_DATABASE}" ]]; then
        echo "PhyloPhlAn database input must be a directory" >&2
        exit 1
    fi

    ln -s "\$PWD/\${NF_DATABASE}" "databases/\${NF_DATABASE_NAME}"

    iqtree_version=\$(iqtree --version 2>&1)
    if ! grep -Eq '(^|[^0-9])3[.]0[.]1([^0-9]|\$)' <<< "\${iqtree_version}"; then
        echo "PhyloPhlAn requires IQ-TREE 3.0.1; configure a combined container or use Conda" >&2
        exit 1
    fi

    echo "Running PhyloPhlAn phylogenomics"

    mkdir -p "\${NF_PREFIX}.phylophlan"

    phylophlan \
        -i genomes \
        -d "\${NF_DATABASE_NAME}" \
        -t a \
        -f "\${NF_CONFIG}" \
        --databases_folder databases \
        --genome_extension .fna \
        --force_nucleotides \
        --diversity high \
        --accurate \
        --nproc "\${NF_CPUS}" \
        --output_folder "\${NF_PREFIX}.phylophlan" \
        --verbose \
        ${args} \
        2>&1 | tee "\${NF_PREFIX}.phylophlan.log"

    mapfile -t tree_files < <(find "\${NF_PREFIX}.phylophlan" -type f -name '*.treefile' -not -path '*/tmp/*')
    if [[ "\${#tree_files[@]}" -ne 1 ]]; then
        echo "Expected one IQ-TREE treefile, found \${#tree_files[@]}" >&2
        exit 1
    fi
    cp "\${tree_files[0]}" "\${NF_PREFIX}.phylophlan.tree.nwk"

    mapfile -t alignments < <(find "\${NF_PREFIX}.phylophlan" -type f -name '*_concatenated.aln' -not -path '*/tmp/*')
    if [[ "\${#alignments[@]}" -ne 1 ]]; then
        echo "Expected one PhyloPhlAn concatenated alignment, found \${#alignments[@]}" >&2
        exit 1
    fi
    cp "\${alignments[0]}" "\${NF_PREFIX}.phylophlan.alignment.fasta"
    test -s "\${NF_PREFIX}.phylophlan.tree.nwk"
    test -s "\${NF_PREFIX}.phylophlan.alignment.fasta"
    rm -rf -- genomes databases
    """

}
