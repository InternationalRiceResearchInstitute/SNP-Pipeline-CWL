#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
-   class: InlineJavascriptRequirement
-   class: StepInputExpressionRequirement
-   class: SubworkflowFeatureRequirement

inputs:
    reference_fa_in: File
    variant_list_in: File[]

outputs:
    vcf_result:
        type: File
        outputSource: genotype_gvcfs/result
steps:
    combine_gvcfs:
        run: ../tools/gatk-combinegvcfs.cwl
        in:
            reference_fa_in: reference_fa_in
            variant_list_in: variant_list_in
            gvcfs_merge_out:
                valueFrom: $(inputs.variant_list_in[0].nameroot.split('_').slice(0)[0])_merged.g.vcf
        out: [result]

    genotype_gvcfs:
        run: ../tools/gatk-genotypegvcfs.cwl
        in:
            reference_fa_in: reference_fa_in
            merged_variant_in: combine_gvcfs/result
            vcf_genotype_out:
                valueFrom: $(inputs.merged_variant_in.nameroot.split('_merged').slice(0)[0]).vcf
        out: [result]


    
doc: |
    GVCFs Merging is a workflow that merges the gvcfs file into one gvcf file and this gvcf file will be converted into vcf file via genotypegvcfs.

    1: MERGE
        GATK v4.0.6.0 --> gatk CombineGVCFs

    2: VCF GENERATION
        GATK v4.0.6.0 --> gatk GenotypeGVCFs
