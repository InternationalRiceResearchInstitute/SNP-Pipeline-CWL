#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement

baseCommand: [java]
arguments: [-Xmx4g, -jar, /usr/local/bin/picard.jar, SortSam]


inputs:
    sam_in:
        type: File
        inputBinding:
            position: 1
            separate: FALSE
            prefix: I=
        doc: <file.sam|file.bam>
    
    bam_sorted_out:
        type: string
        default: sorted.bam
        inputBinding:
            position: 2
            separate: FALSE
            prefix: O=
        doc: Filename for the output sorted bam file
    
    sort_order:
        type: string
        default: coordinate
        inputBinding:
            position: 3
            separate: FALSE
            prefix: SO=
    
    validation_stringency:
        type: string
        default: LENIENT
        inputBinding:
            position: 4
            separate: FALSE
            prefix: VALIDATION_STRINGENCY=

outputs:
    result:
        type: File
        outputBinding:
            glob: $(inputs.bam_sorted_out)
        
doc: |
    USAGE: SortSam [options]

    Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#SortSam

    This tool sorts the input SAM or BAM file by coordinate, queryname (QNAME), or some other property of the SAM record.
    The SortOrder of a SAM/BAM file is found in the SAM file header tag @HD in the field labeled SO.  
    For a coordinate sorted SAM/BAM file, read alignments are sorted first by the reference sequence name (RNAME) field
    using the reference sequence dictionary (@SQ tag).  Alignments within these subgroups are secondarily sorted using the
    left-most mapping position of the read (POS).  Subsequent to this sorting scheme, alignments are listed arbitrarily.

    For queryname-sorted alignments, all alignments are grouped using the queryname field but the alignments are not
    necessarily sorted within these groups.  Reads having the same queryname are derived from the same template.

    Usage example:

    java -jar picard.jar SortSam \
    I=input.bam \
    O=sorted.bam \
    SORT_ORDER=coordinate

    Version: 2.18.9-SNAPSHOT
