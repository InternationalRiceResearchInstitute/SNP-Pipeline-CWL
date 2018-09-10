#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
-   class: InlineJavascriptRequirement

baseCommand: [gatk, GenotypeGVCFs]

inputs:
    merged_variant_in:
        type: File
        inputBinding:
            position: 1
            prefix: --variant
        doc: <file.g.vcf>
    
    vcf_genotype_out:
        type: string
        default: vcf_genotype_out.vcf
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
            glob: $(inputs.vcf_genotype_out)
        
doc: |
    USAGE: GenotypeGVCFs [arguments]

    Perform joint genotyping on a single-sample GVCF from HaplotypeCaller or a multi-sample GVCF from CombineGVCFs or
    GenomicsDBImport
    Version:4.0.6.0


    Required Arguments:

    --output,-O:File              File to which variants should be written  Required. 

    --reference,-R:String         Reference sequence file  Required. 

    --variant,-V:String           A VCF file containing variants  Required. 

