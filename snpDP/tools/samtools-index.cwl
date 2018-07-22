#!/usr/bin/env cwl-runnner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement
-   class: InitialWorkDirRequirement
    listing:
    -   entry: $(inputs.bam_in)
        entryname: $(inputs.bam_in.path.split('/').slice(-1)[0])

baseCommand: [samtools, index]

inputs: 
    bam_in:
        type: File
        inputBinding:
            position: 1

outputs:
    result: 
        type: File
        secondaryFiles:
            -   .bai
        outputBinding:
            glob: $(inputs.bam_in.path.split('/').slice(-1)[0])
    