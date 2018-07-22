#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement

baseCommand: [java]
arguments: [-Xmx4g, -jar, /usr/local/bin/picard.jar, FixMateInformation]


inputs:
    sorted_in:
        type: File
        inputBinding:
            position: 1
            separate: FALSE
            prefix: I=
        doc: <file.sam|file.bam>
    
    fixmateinfo_out:
        type: string
        default: fxmt.bam
        inputBinding:
            position: 2
            separate: FALSE
            prefix: O=
        doc: Filename for the output fix mate bam file
    
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
            glob: $(inputs.fixmateinfo_out)
        
doc: |
    USAGE: FixMateInformation [options]

    Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#FixMateInformation

    Verify mate-pair information between mates and fix if needed.This tool ensures that all mate-pair information is in sync
    between each read and its mate pair.  If no OUTPUT file is supplied then the output is written to a temporary file and
    then copied over the INPUT file (with the original placed in a .old file.)  Reads marked with the secondary alignment
    flag are written to the output file unchanged. However supplementary</b> reads are corrected so that they point to the
    primary, non-supplemental mate record.


    Usage example


    java -jar picard.jar FixMateInformation \
    I=input.bam \
    O=fixed_mate.bam \
    ADD_MATE_CIGAR=true


    Caveats

    The program should run with fairly limited memory unless there are many mate pairs that are missing or far apart from
    each other in the file, as it keeps track of the unmatched mates.
    Version: 2.18.9-SNAPSHOT
