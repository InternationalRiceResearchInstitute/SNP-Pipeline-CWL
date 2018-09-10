#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement
-   class: InitialWorkDirRequirement
    listing:
    -   entry: $(inputs.reference_fa_in)
        entryname: $(inputs.reference_fa_in.path.split('/').slice(-1)[0])

baseCommand: [java]
arguments: [-Xmx4g, -jar, /usr/local/bin/picard.jar, CreateSequenceDictionary]


inputs:
    reference_fa_in:
        type: File
        secondaryFiles:
            -    .fai
        inputBinding:
            position: 1
            prefix: R=
        doc: <file.fasta|file.fa|file.fa.gz>
    
    dict_out:
        type: string
        default: ref.dict
        inputBinding:
            position: 2
            prefix: O=
        doc: Filename for the output dictionary file
    
outputs:
    dict_result:
        type: File
        secondaryFiles:
            -   .fai
            -   '^.dict'
        outputBinding:
            glob: $(inputs.reference_fa_in.basename)
        
doc: |
    USAGE: CreateSequenceDictionary [options]

    Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary

    Creates a sequence dictionary for a reference sequence.  This tool creates a sequence dictionary file (with ".dict"
    extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools.
    The output file contains a header but no SAMRecords, and the header contains only sequence records.

    The reference sequence can be gzipped (both .fasta and .fasta.gz are supported).
    Usage example:

    java -jar picard.jar CreateSequenceDictionary \ 
    R=reference.fasta \ 
    O=reference.dict

    Version: 2.18.9-SNAPSHOT
