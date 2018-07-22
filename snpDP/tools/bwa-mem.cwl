#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement
-   class: InitialWorkDirRequirement
    listing:
    -   entry: $(inputs.sam_out)
        entryname: $(inputs.sam_out)

baseCommand: [bwa, mem, -M]

inputs:
    bwa_threads:
        type: int
        default: 8
        inputBinding:
            position: 1
            prefix: -t
        doc: Number of threads to be used for bwa mem tool

    reference_fa_in:
        type: File
        secondaryFiles:
            -   .amb
            -   .ann
            -   .bwt
            -   .pac
            -   .sa
        inputBinding:
            position: 2
        doc: <file.fasta|file.fa|file.fa.gz>

    reads1_in:
        type: File
        inputBinding:
            position: 3
        doc: <file.fastq|file.fq|file.fq.gz>
    
    reads2_in:
        type: File
        inputBinding:
            position: 4
        doc: <file.fastq|file.fq|file.fq.gz>
    
    sam_out:
        type: string
        default: sam_out.sam

stdout: $(inputs.sam_out)

outputs:
    result:
        type: stdout
        
doc: |
    Usage:   bwa index [options] <in.fasta>

    Options:    -a STR    BWT construction algorithm: bwtsw or is [auto]
                -p STR    prefix of the index [same as fasta name]
                -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
                -6        index files named as <in.fasta>.64.* instead of <in.fasta>.*

    Warning:    `-a bwtsw' does not work for short genomes, while `-a is' and
                `-a div' do not work not for long genomes.