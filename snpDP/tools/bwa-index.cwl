#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement
-   class: InitialWorkDirRequirement
    listing:
    -   entry: $(inputs.reference_fa_in)
        entryname: $(inputs.reference_fa_in.path.split('/').slice(-1)[0])

baseCommand: [bwa, index]

inputs:
    reference_fa_in:
        type: File
        secondaryFiles:
            -    .fai
            -   '^.dict'
        inputBinding:
            position: 1
        doc: <file.fasta|file.fa|file.fa.gz>
        
outputs:
    bwa_indices_result:
        type: File
        secondaryFiles:
            -   .fai
            -   '^.dict'
            -   .amb
            -   .ann
            -   .bwt
            -   .pac
            -   .sa
        outputBinding:
            glob: $(inputs.reference_fa_in.path.split('/').slice(-1)[0])
        
doc: |
    Usage:   bwa index [options] <in.fasta>

    Options:    -a STR    BWT construction algorithm: bwtsw or is [auto]
                -p STR    prefix of the index [same as fasta name]
                -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
                -6        index files named as <in.fasta>.64.* instead of <in.fasta>.*

    Warning:    `-a bwtsw' does not work for short genomes, while `-a is' and
                `-a div' do not work not for long genomes.