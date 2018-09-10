#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement

baseCommand: [gatk, CombineGVCFs]

inputs:
    variant_list_in:
        type:
            type: array
            items: File
            inputBinding:
                position: 1
                prefix: --variant
        doc: <file.g.vcf>
    
    gvcfs_merge_out:
        type: string
        default: gvcfs_merge_out.g.vcf
        inputBinding:
            position: 2
            prefix: --output
        doc: Filename for the variant calling output gvcf file
    
    reference_fa_in:
        type: File
        inputBinding:
            position: 3
            prefix: --reference

outputs:
    result:
        type: File
        outputBinding:
            glob: $(inputs.gvcfs_merge_out)
        
doc: |
    USAGE: CombineGVCFs [arguments]

    Merges one or more HaplotypeCaller GVCF files into a single GVCF with appropriate annotations
    Version:4.0.6.0


    Required Arguments:

    --output,-O:File              The combined GVCF output file  Required. 

    --reference,-R:String         Reference sequence file  Required. 

    --variant,-V:String           One or more VCF files containing variants  This argument must be specified at least once.
                                Required. 
