#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement

baseCommand: [samtools, merge]

inputs:
    bam_merge_out:
        type: string
        default: bam_merged_out.bam
        inputBinding:
            position: 1
        doc: Output filename for the merged bam file
    
    bam_list_in:
        type: File[]
        inputBinding:
            position: 2
        doc: List of bam files 

outputs:
    result:
        type: File
        outputBinding:
            glob: $(inputs.bam_merge_out)
        
doc: |
    Usage: samtools merge [-nurlf] [-h inh.sam] [-b <bamlist.fofn>] <out.bam> <in1.bam> [<in2.bam> ... <inN.bam>]
