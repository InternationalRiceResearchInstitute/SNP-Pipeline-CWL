#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement

baseCommand: [java]
arguments: [-Xmx4g, -jar, /usr/local/bin/picard.jar, MarkDuplicates]


inputs:
    fixmate_in:
        type: File
        inputBinding:
            position: 1
            separate: FALSE
            prefix: I=
        doc: <file.bam>
    
    markdup_out:
        type: string
        default: mkdup.bam
        inputBinding:
            position: 2
            separate: FALSE
            prefix: O=
        doc: Filename for the output mark duplicate bam file
    
    validation_stringency:
        type: string
        default: LENIENT
        inputBinding:
            position: 3
            separate: FALSE
            prefix: VALIDATION_STRINGENCY=
    
    metrics_out:
        type: string
        default: out.metrics
        inputBinding:
            position: 4
            separate: FALSE
            prefix: METRICS_FILE=

outputs:
    result:
        type: File
        outputBinding:
            glob: $(inputs.markdup_out)
        
doc: |
    USAGE: MarkDuplicates [options]

    Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

    Identifies duplicate reads.  
    This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a
    single fragment of DNA.  Duplicates can arise during sample preparation e.g. library construction using PCR.  See also
    EstimateLibraryComplexity (https://broadinstitute.github.io/picard/command-line-overview.html#EstimateLibraryComplexity)
    for additional notes on PCR duplication artifacts.  Duplicate reads can also result from a single amplification cluster,
    incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.  These duplication
    artifacts are referred to as optical duplicates.

    The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM
    file.  An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes.  After duplicate
    reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the
    sums of their base-quality scores (default method).

    The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for
    each read.  Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024. 
    If you are not familiar with this type of annotation, please see the following blog post
    (https://www.broadinstitute.org/gatk/blog?id=7019) for additional information.

    Although the bitwise flag annotation indicates whether a read was marked as a duplicate, it does not identify the type
    of duplicate.  To do this, a new tag called the duplicate type (DT) tag was recently added as an optional output in  the
    'optional field' section of a SAM/BAM file.  Invoking the TAGGING_POLICY option, you can instruct the program to mark
    all the duplicates (All), only the optical duplicates (OpticalOnly), or no duplicates (DontTag).  The records within the
    output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), as either
    library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ).  This tool uses the
    READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the primary methods to identify and differentiate
    duplicate types.  Set READ_NAME_REGEX to null to skip optical duplicate detection, e.g. for RNA-seq or other data where
    duplicate sets are extremely large and estimating library complexity is not an aim.  Note that without optical duplicate
    counts, library size estimation will be inaccurate.

    MarkDuplicates also produces a metrics file indicating the numbers of duplicates for both single- and paired-end reads.

    The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different.  When
    the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked
    as duplicates.  However, when the input is query-sorted (actually query-grouped), then unmapped mates and
    secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.

    If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.

    Usage example:

    java -jar picard.jar MarkDuplicates \
    I=input.bam \
    O=marked_duplicates.bam \
    M=marked_dup_metrics.txt
    Please see MarkDuplicates (http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics) for
    detailed explanations of the output metrics.

    Version: 2.18.9-SNAPSHOT

