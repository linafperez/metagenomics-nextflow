/** Input validation and metadata serialization shared by both input modes. */
class InputMetadata {
    static List<Map> table(List records, List<String> required, String groupColumn = '') {
        if (!records || records.size() < 2) {
            throw new IllegalArgumentException('Input table requires a header and at least one sample')
        }
        def header = records[0].collect { it == null ? '' : it.toString().trim() }
        header[0] = header[0].replaceFirst(/^\uFEFF/, '')
        if (header.any { !it } || header.toSet().size() != header.size()) {
            throw new IllegalArgumentException('Input table has empty or duplicate column names')
        }
        def missing = (required + (groupColumn ? [groupColumn] : [])).findAll { !header.contains(it) }
        if (missing) {
            throw new IllegalArgumentException("Input table is missing columns: ${missing.join(', ')}")
        }
        if (groupColumn && groupColumn in required) {
            throw new IllegalArgumentException('An identity/read column cannot also be the group column')
        }
        records.drop(1).withIndex().collect { values, index ->
            if (values.size() != header.size()) {
                throw new IllegalArgumentException("Input line ${index + 2} has ${values.size()} fields; expected ${header.size()}")
            }
            def row = [:]
            header.eachWithIndex { key, i ->
                def value = values[i] == null ? '' : values[i].toString().trim()
                if (value.find(/[\t\r\n]/)) {
                    throw new IllegalArgumentException("Input line ${index + 2}, column ${key} contains a tab or line break")
                }
                row[key] = value
            }
            row
        }
    }

    static void validateSamples(List<Map> samples, boolean condition) {
        def ids = [] as Set
        samples.each { meta ->
            if (!(meta.id ==~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/) || !ids.add(meta.id.toLowerCase(Locale.ROOT))) {
                throw new IllegalArgumentException("Invalid or duplicate sample ID: ${meta.id}")
            }
        }
        if (condition) {
            // Validate groups before any download or preprocessing task is submitted.
            SpadesCoassembly.conditionInputs(samples.collect { [it, ['mate1', 'mate2']] })
        }
    }

    static Map sample(Map row, String id, String groupColumn, String project = '') {
        def groupKey = groupColumn ?: (row.containsKey('group') ? 'group' : '')
        def group = groupKey ? row[groupKey] : ''
        if (groupKey && !group) {
            throw new IllegalArgumentException("Sample ${id} has an empty ${groupKey} group")
        }
        // Preserve every input column, including study-specific metadata.
        [id: id, single_end: false, biosample_accession: project ? id : '',
         identity_source: project ? 'BioSample' : 'samplesheet', group: group,
         project_accession: project, run_accessions: [], metadata: new LinkedHashMap(row)]
    }

    static List<Map> sraSamples(List records, String groupColumn, String project, boolean condition) {
        def rows = table(records, ['biosample_accession'], groupColumn)
        def samples = rows.withIndex().collect { row, i ->
            def accession = row.biosample_accession.toUpperCase(Locale.ROOT)
            if (!(accession ==~ /^SAM(?:N|EA|D)[0-9]+$/)) {
                throw new IllegalArgumentException("Invalid BioSample accession: ${accession}")
            }
            sample(row, accession, groupColumn, project)
        }
        validateSamples(samples, condition)
        samples
    }

    static List<Map> localSamples(List records, String groupColumn, boolean condition) {
        def samples = table(records, ['sample', 'fastq_1', 'fastq_2'], groupColumn)
            .collect { row -> sample(row, row.sample, groupColumn) }
        validateSamples(samples, condition)
        samples
    }

    static String metadataTsv(List<Map> samples) {
        def fields = ['sample_id', 'biosample_accession', 'group', 'run_accessions', 'project_accession', 'sample_order']
        def extra = samples.collectMany { it.metadata.keySet().toList() }.unique().sort()
            .findAll { !(it in fields) }
        def header = fields + extra
        def rows = samples.withIndex().collect { meta, index ->
            def row = meta.metadata + [sample_id: meta.id, biosample_accession: meta.biosample_accession,
                group: meta.group, run_accessions: meta.run_accessions.join(';'),
                project_accession: meta.project_accession, sample_order: index + 1]
            header.collect { row[it]?.toString() ?: '' }.join('\t')
        }
        ([header.join('\t')] + rows).join('\n') + '\n'
    }
}
