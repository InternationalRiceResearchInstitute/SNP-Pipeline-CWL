#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement
-   class: InitialWorkDirRequirement
    listing:
    -   entry: $(inputs.reference_fa_in)
        entryname: $(inputs.reference_fa_in.path.split('/').slice(-1)[0])

baseCommand: [samtools, faidx]

inputs:
    reference_fa_in:
        type: File
        inputBinding:
            position: 1
        doc: <file.fasta|file.fa|file.fa.gz>
        
outputs:
    fai_result:
        type: File
        secondaryFiles:
            -   .fai
        outputBinding:
            glob: $(inputs.reference_fa_in.path.split('/').slice(-1)[0])
        
doc: |
    samtools-faidx.cwl is developed for CWL consortium
    Usage:   samtools faidx <file.fa|file.fa.gz> [<reg> [...]]
