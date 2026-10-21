# Streaming validation of uncompressed, four-line paired FASTQ files.
# Usage: awk -v mate2=R2.fastq -f validate_fastq_pair.awk R1.fastq
function fail(message) {
    print "Invalid paired FASTQ: " message > "/dev/stderr"
    failed = 1
    exit 1
}
function read_id(header, mate, fields, count, name) {
    if (substr(header, 1, 1) != "@") fail("missing @ header")
    count = split(substr(header, 2), fields, /[ \t]+/)
    name = fields[1]
    if (!length(name)) fail("empty read ID")
    if (name ~ /\/[12]$/) {
        if (substr(name, length(name), 1) != mate) fail("wrong mate suffix")
        sub(/\/[12]$/, "", name)
    }
    if (count > 1 && fields[2] ~ /^[12]:/ && substr(fields[2], 1, 1) != mate) fail("wrong CASAVA mate")
    return name
}
{
    sub(/\r$/, "")
    if ((getline second < mate2) != 1) fail("mate 2 is truncated")
    sub(/\r$/, "", second)
    position = (NR - 1) % 4
    if (position == 0 && read_id($0, 1) != read_id(second, 2)) fail("different read IDs at record " int((NR + 3) / 4))
    if (position == 1) {
        length1 = length($0)
        length2 = length(second)
        if (!length1 || !length2) fail("empty sequence")
    }
    if (position == 2 && (substr($0, 1, 1) != "+" || substr(second, 1, 1) != "+")) fail("missing + separator")
    if (position == 3 && (length($0) != length1 || length(second) != length2)) fail("sequence/quality length mismatch")
}
END {
    if (failed) exit 1
    if (!NR || NR % 4) fail("empty or incomplete records")
    if ((getline second < mate2) != 0) fail("mate 2 has extra records or an I/O error")
    close(mate2)
    print "Validated " NR / 4 " paired records"
}
