/** Resolve an explicitly selected cohort against public NCBI RunInfo metadata. */
class SraMetadata {
    static String field(Map row, String key) {
        def aliases = [
            run: ['Run', 'run_accession'], biosample: ['BioSample', 'biosample_accession'],
            project: ['BioProject', 'project_accession'], experiment: ['Experiment', 'experiment_accession'],
            layout: ['LibraryLayout', 'layout'], strategy: ['LibraryStrategy', 'strategy'],
            source: ['LibrarySource', 'source'], platform: ['Platform', 'platform'],
            download: ['download_path', 'DownloadPath'], consent: ['Consent', 'consent'],
            access: ['public_access', 'Access', 'access'], dbgap: ['dbgap_study_accession', 'dbGaP'],
            spots: ['spots'], mates: ['spots_with_mates', 'spotsWithMates'], bases: ['bases']
        ]
        def normalized = row.collectEntries { k, v -> [(k.toLowerCase(Locale.ROOT).replaceAll(/[^a-z0-9]/, '')): v] }
        def value = aliases[key].collect { normalized[it.toLowerCase(Locale.ROOT).replaceAll(/[^a-z0-9]/, '')] }
            .find { it != null && it.toString().trim() }
        value?.toString()?.trim() ?: ''
    }

    static Map resolve(List<Map> samples, List records, String project, String platforms) {
        def rows = InputMetadata.table(records, [])
        def headers = records[0].collect { it.toString().toLowerCase(Locale.ROOT).replaceAll(/[^a-z0-9]/, '') }
        if (headers.toSet().size() != headers.size()) {
            throw new IllegalArgumentException('RunInfo has duplicate column names after normalization')
        }
        def selected = samples.collectEntries { [(it.id): it] }
        def allowed = platforms.tokenize(',').collect { it.trim().toUpperCase(Locale.ROOT) } as Set
        def seen = [] as Set
        def eligible = [:].withDefault { [] }
        def exclusions = []
        rows.each { row ->
            def biosample = field(row, 'biosample').toUpperCase(Locale.ROOT)
            if (selected.containsKey(biosample)) {
                def run = field(row, 'run').toUpperCase(Locale.ROOT)
                if (!seen.add(run)) {
                    throw new IllegalArgumentException("Duplicate or contradictory selected RunInfo accession: ${run}")
                }
                def reasons = []
                if (!(run ==~ /^[SED]RR[0-9]+$/)) reasons << 'invalid_run_accession'
                if (!field(row, 'project').toUpperCase(Locale.ROOT).findAll(/PRJ(?:NA|EB|DB)[0-9]+/).contains(project)) reasons << 'bioproject_accession_mismatch'
                if (field(row, 'layout').toUpperCase(Locale.ROOT) != 'PAIRED') reasons << 'library_layout_not_paired'
                if (field(row, 'strategy').toUpperCase(Locale.ROOT) != 'WGS') reasons << 'library_strategy_not_wgs'
                if (field(row, 'source').toUpperCase(Locale.ROOT) != 'METAGENOMIC') reasons << 'library_source_not_metagenomic'
                if (!allowed.contains(field(row, 'platform').toUpperCase(Locale.ROOT))) reasons << 'unsupported_platform'
                def consent = field(row, 'consent').toLowerCase(Locale.ROOT)
                def publicFlag = field(row, 'access').toLowerCase(Locale.ROOT)
                def download = field(row, 'download')
                def denied = ['false', 'no', '0', 'controlled', 'restricted', 'protected']
                if (publicFlag in denied || consent in denied || field(row, 'dbgap') || !download ||
                    (publicFlag && !(publicFlag in ['true', 'yes', '1', 'public', 'open']))) {
                    reasons << 'non_public_or_unknown_access'
                }
                def counts = [:]
                ['spots', 'mates', 'bases'].each { key ->
                    def raw = field(row, key)
                    if (raw) {
                        if (!(raw ==~ /^[+-]?[0-9]+$/)) {
                            reasons << "invalid_${key}"
                        } else {
                            counts[key] = new BigInteger(raw)
                            if (key != 'mates' && counts[key] <= 0) reasons << "non_positive_${key}"
                            if (key == 'mates' && counts[key] < 0) reasons << 'negative_spots_with_mates'
                        }
                    }
                }
                if (counts.containsKey('spots') && counts.containsKey('mates') && counts.mates > counts.spots) {
                    reasons << 'spots_with_mates_exceeds_spots'
                }
                // Fewer paired spots than spots is advisory; validate actual FASTQ pairs in the acquisition task.
                if (reasons) {
                    exclusions << [biosample, run, reasons.join(';')]
                } else {
                    eligible[biosample] << run
                }
            }
        }
        def resolved = samples.collect { meta ->
            def runs = eligible[meta.id].sort()
            if (!runs) {
                throw new IllegalArgumentException("BioSample ${meta.id} has no public, paired WGS METAGENOMIC runs in ${project} on ${platforms}")
            }
            meta + [run_accessions: runs]
        }
        [samples: resolved, exclusions: exclusions]
    }

    static String exclusionsTsv(Map resolution) {
        (['biosample_accession\trun_accession\texclusion_reason'] +
            resolution.exclusions.collect { it.join('\t') }).join('\n') + '\n'
    }
}
