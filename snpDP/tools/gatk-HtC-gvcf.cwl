#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement

baseCommand: [gatk, HaplotypeCaller]

inputs:
    addrep_in:
        type: File
        secondaryFiles:
            -   '^.bai'
        inputBinding:
            position: 1
            prefix: --input
        doc: <file.bam and file.bai>
    
    gvcf_out:
        type: string
        default: gvcf_out.g.vcf
        inputBinding:
            position: 2
            prefix: --output
        doc: Filename for the variant calling output gvcf file
    
    reference_fa_in:
        type: File
        inputBinding:
            position: 3
            prefix: --reference

    emit_ref_conf:
        type: string
        default: GVCF
        inputBinding:
            position: 4
            prefix: --emit-ref-confidence

outputs:
    result:
        type: File
        outputBinding:
            glob: $(inputs.gvcf_out)
        
doc: |
    USAGE: HaplotypeCaller [arguments]

    Call germline SNPs and indels via local re-assembly of haplotypes
    Version:4.0.6.0


    Required Arguments:

    --input,-I:String             BAM/SAM/CRAM file containing reads  This argument must be specified at least once.
                                Required. 

    --output,-O:String            File to which variants should be written  Required. 

    --reference,-R:String         Reference sequence file  Required.

