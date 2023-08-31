import pysam


# pay attention which input is specified
# 'fillTheGaps_test.sorted.bam' is for the testing
# 'Aste.sorted.bam' is for the data processing
samfile = pysam.AlignmentFile("Aste.sorted.bam", "rb")
acceptedReads = []
gapCoord = [14928794, 14931978]  # these are certain gap borders (#2.5); in the future there should be a list of gaps designated manually
clippedReads = open('fillTheGaps_2.5_output.fa', 'w', encoding='utf8')

# specify chromosome name in parentheses
for read in samfile.fetch('X'):
    if read.is_mapped and not read.is_supplementary:
        refStart = read.reference_start
        qStart = read.query_alignment_start
        coord_chrom = refStart - qStart  # to define the beginning of a read in the coordinate system of a chromosome
        qLen = read.query_length

        # the beginning of a read should be before the end of the gap &
        # the end of a read should be after the beginning of the gap
        if coord_chrom < gapCoord[1] and qLen + coord_chrom > gapCoord[0]:
            acceptedReads.append((read.query_name,
                                  '+' if read.is_forward else '-',
                                  refStart,
                                  read.reference_end,
                                  qStart,
                                  read.query_alignment_end,
                                  read.query_sequence)
                                 )

# output in FASTA format
for read in acceptedReads:
    header = '>' + '\t'.join(map(str, read[0:6])) + '\n'
    coord_chrom = read[2] - read[4]
    coord_read = read[4] - read[2]
    read_sequence = read[6]

    # if the end of a read is after the end of the gap, it must be sliced on the right
    if len(read_sequence) + coord_chrom > gapCoord[1]:
        read_sequence = read_sequence[:gapCoord[1] + coord_read + 1]

    # if the beginning of a read is before the beginning of the gap, it must be sliced on the left
    if coord_chrom < gapCoord[0]:
        clippedReads.write(header + read_sequence[gapCoord[0] + coord_read:] + '\n')
    else:
        clippedReads.write(header + read_sequence + '\n')

    # the situation when a read is entirely inside the gap should not happen —
    # it means that this read is unmapped, which contradicts with the 'read.is_mapped' condition

clippedReads.close()
