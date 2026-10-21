# Filter ASCII FASTA records; write FASTA to stdout and metrics to stats_file.
# The caller publishes stdout only after successful validation.
BEGIN {
    RS = "\r\n|\r|\n"
    if (minimum_length !~ /^[0-9]+$/ || minimum_length + 0 < 1) {
        fail("Minimum length must be a positive integer")
    }
    if (stats_file == "") fail("stats_file is required")
    total_records = total_bases = kept_records = kept_bases = 0
}

function fail(message) {
    print "FASTA filter: " message > "/dev/stderr"
    failed = 1
    exit 1
}

function emit_record(    size, offset) {
    if (header == "") return
    size = length(sequence)
    total_records++
    total_bases += size
    if (size < minimum_length) return
    kept_records++
    kept_bases += size
    print header
    for (offset = 1; offset <= size; offset += 80) {
        print substr(sequence, offset, 80)
    }
}

{
    line = $0
    sub(/^[[:space:]]+/, "", line)
    sub(/[[:space:]]+$/, "", line)
    if (line == "") next
    if (substr(line, 1, 1) == ">") {
        emit_record()
        header = line
        sequence = ""
    } else {
        if (header == "") fail("Sequence found before FASTA header")
        sequence = sequence line
    }
}

END {
    if (failed) exit 1
    emit_record()
    if (kept_records == 0) fail("No contigs meet the " minimum_length " bp threshold")
    print "metric\tvalue" > stats_file
    printf "minimum_length\t%.0f\n", minimum_length > stats_file
    printf "input_contigs\t%.0f\n", total_records > stats_file
    printf "input_bases\t%.0f\n", total_bases > stats_file
    printf "retained_contigs\t%.0f\n", kept_records > stats_file
    printf "retained_bases\t%.0f\n", kept_bases > stats_file
    if (close(stats_file) != 0) fail("Cannot write statistics")
}
