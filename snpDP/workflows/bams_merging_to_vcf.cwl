#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
-   class: InlineJavascriptRequirement
-   class: StepInputExpressionRequirement
-   class: SubworkflowFeatureRequirement

inputs:
    reference_fa_in: File
    bam_list_in: File[]

outputs:
    vcf_result:
        type: File
        outputSource: gatk_HtC_vcf/result
steps:
    samtools_merge:
        run: ../tools/samtools-merge.cwl
        in:
            bam_list_in: bam_list_in
            bam_merge_out:
                valueFrom: $(inputs.bam_list_in[0].nameroot.split('_').slice(0)[0])_merged.bam
        out: [result]
    
    samtools_index:
        run: ../tools/samtools-index.cwl
        in:
            bam_in: samtools_merge/result
        out: [result]

    gatk_HtC_vcf:
        run: ../tools/gatk-HtC-vcf.cwl
        in:
            reference_fa_in: reference_fa_in
            bam_in: samtools_index/result
            vcf_out:
                valueFrom: $(inputs.bam_in.nameroot.split('_merged').slice(0)[0]).vcf
        out: [result]

    
doc: |
    GVCFs Merging is a workflow that merges the gvcfs file into one gvcf file and this gvcf file will be converted into vcf file via genotypegvcfs.

    1: MERGE
        GATK v4.0.6.0 --> gatk CombineGVCFs

    2: VCF GENERATION
        GATK v4.0.6.0 --> gatk GenotypeGVCFs
