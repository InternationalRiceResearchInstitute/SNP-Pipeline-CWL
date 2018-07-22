#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement

baseCommand: [gatk, HaplotypeCaller]

inputs:
    reference_fa_in:
        type: File
        secondaryFiles:
            -   .fai
            -   '^.dict'
        inputBinding:
            position: 1
            prefix: -R
    
    bam_in:
        type: File
        secondaryFiles:
            - .bai
        inputBinding:
            position: 2
            prefix: -I
    
    vcf_out:
        type: string
        default: gatk_HtC_vcf_out.vcf
        inputBinding:
            position: 3
            prefix: -O

outputs:
    result:
        type: File
        outputBinding:
            glob: $(inputs.vcf_out)