/** Condition-only read validation and deterministic filesystem identifiers. */
class SpadesCoassembly {
    // Keep this rule identical to group_id() in bin/spades_coassembly.py.
    static String groupId(String group) {
        'group_' + group.replaceAll(/[^A-Za-z0-9._-]+/, '_').take(80)
    }

    static List conditionInputs(List records) {
        if (!records) {
            throw new IllegalArgumentException('SPAdes condition coassembly requires at least one paired-end sample')
        }
        def seenSamples = [] as Set
        def safeGroups = [:]
        records.each { record ->
            if (record.size() != 2 || !record[0]['id'] || record[1].size() != 2) {
                throw new IllegalArgumentException('SPAdes coassembly expects tuples of [meta, [read_1, read_2]]')
            }
            def id = record[0]['id'].toString()
            if (!seenSamples.add(id)) {
                throw new IllegalArgumentException("SPAdes condition coassembly has duplicate/ambiguous sample ID: ${id}")
            }
            def group = record[0]['group']?.toString()
            if (!group || !group.trim() || group.find(/[\t\r\n]/)) {
                throw new IllegalArgumentException("SPAdes condition coassembly requires a non-empty group for sample ${id}")
            }
            def key = groupId(group).toLowerCase(java.util.Locale.ROOT)
            if (safeGroups.containsKey(key) && safeGroups[key] != group) {
                throw new IllegalArgumentException("SPAdes group filename collision: '${safeGroups[key]}' and '${group}'")
            }
            safeGroups[key] = group
        }
        records.groupBy { record -> record[0]['group'].toString() }
            .entrySet().toList().sort { left, right -> left.key <=> right.key }
            .collect { entry ->
                def ordered = entry.value.toList().sort { left, right ->
                    left[0]['id'].toString() <=> right[0]['id'].toString()
                }
                def safe = groupId(entry.key)
                def meta = [
                    id: "spades_coassembly_${safe}".toString(),
                    assembler: 'spades', branch: 'spades',
                    coassembly_mode: 'condition', group: entry.key, group_id: safe,
                    sample_ids: ordered.collect { record -> record[0]['id'].toString() }
                ]
                [meta, ordered.collectMany { record -> record[1].toList() }]
            }
    }
}
