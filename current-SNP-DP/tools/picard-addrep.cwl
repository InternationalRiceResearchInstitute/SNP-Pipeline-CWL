#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement

baseCommand: [java]
arguments: [-Xmx4g, -jar, /usr/local/bin/picard.jar, AddOrReplaceReadGroups]


inputs:
    markdup_in:
        type: File
        inputBinding:
            position: 1
            separate: FALSE
            prefix: I=
        doc: <file.bam>
    
    addrep_out:
        type: string
        default: addrep.bam
        inputBinding:
            position: 2
            separate: FALSE
            prefix: O=
        doc: Filename for the output addrep bam file
    
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
    
    create_index:
        type: string
        default: 'TRUE'
        inputBinding:
            position: 5
            separate: FALSE
            prefix: CREATE_INDEX=

    rgid: 
        type: string
        default: addrep_id
        inputBinding: 
            position: 6
            separate: FALSE
            prefix: RGID=

    rglb:
        type: string
        default: addrep_lb
        inputBinding: 
            position: 6
            separate: FALSE
            prefix: RGLB=

    rgpl:
        type: string
        default: ILLUMINA
        inputBinding: 
            position: 6
            separate: FALSE
            prefix: RGPL=
    rgpu:
        type: string
        default: addrep_pu
        inputBinding: 
            position: 6
            separate: FALSE
            prefix: RGPU=
    rgsm: 
        type: string
        default: addrep_sm
        inputBinding: 
            position: 6
            separate: FALSE
            prefix: RGSM=

outputs:
    result:
        type: File
        secondaryFiles:
            -   '^.bai'
        outputBinding:
            glob: $(inputs.addrep_out)
        
doc: |
    USAGE: AddOrReplaceReadGroups [options]

    Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups

    Assigns all the reads in a file to a single new read-group.

    This tool accepts INPUT BAM and SAM files or URLs from the Global Alliance for Genomics and Health (GA4GH)
    (http://ga4gh.org/#/documentation).

    Usage example:

    java -jar picard.jar AddOrReplaceReadGroups \
    I=input.bam \
    O=output.bam \
    RGID=4 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=20


    Caveats

    The value of the tags must adhere (according to the SAM-spec (https://samtools.github.io/hts-specs/SAMv1.pdf)) with the
    regex '^[ -~]+$'</code> (one or more characters from the ASCII range 32 through 126). In particular <Space> is the only
    non-printing character allowed.

    The program enables only the wholesale assignment of all the reads in the INPUT to a single read-group. If your file
    already has reads assigned to multiple read-groups, the original RG value will be lost. 

    For more information about read-groups, see the GATK Dictionary entry.
    (https://www.broadinstitute.org/gatk/guide/article?id=6472)
    Version: 2.18.9-SNAPSHOT


