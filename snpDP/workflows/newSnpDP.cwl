#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
-   class: SubworkflowFeatureRequirement
-   class: ScatterFeatureRequirement

inputs:
    reference_fa_in: File
    reads1_in: File[]
    reads2_in: File[]

outputs:
    new_snpDP_result:
        type: File
        outputSource: gvcfs_merging/vcf_result

steps:
    indices_generation:
        run: indices_generation.cwl
        in:
            reference_fa_in: reference_fa_in
        out: [indices_result]

    gvcfs_generation:
        scatter: [reads1_in, reads2_in]
        scatterMethod: dotproduct
        run: gvcfs_generation.cwl
        in: 
            reference_fa_in: indices_generation/indices_result
            reads1_in: reads1_in
            reads2_in: reads2_in
        out: [result]
    
    gvcfs_merging:
        run: gvcfs_merging.cwl
        in:
            reference_fa_in: indices_generation/indices_result
            variant_list_in: gvcfs_generation/result
        out: [vcf_result]

doc: |
    New SNP-Discovery-Pipeline

    The New SNP discovery pipeline like the current SNP discovery pipeline is used to detect variants including SNPs and indels (insertion/deletions) from next-generation sequencing (NGS) reads. Instead having the step MERGE BAM then VARIANT CALLING smaller BAM files fed VARIANT CALLING then MERGE GVCF output files to CREATE VCF. The pipeline includes the following:

    Alignment: The sequence reads are aligned to the reference genome using BWA (Burrows-Wheeler Aligner) – bwa 0.7.17 http://bio-bwa.sourceforge.net/

    BAM Processing: A series of intermediate steps to process and prepare the BAM file for variant calling. Picard Tools is used for BAM processing – Picard Tools 2.18.9 http://broadinstitute.github.io/picard/

    Variant calling: Variants are called using the GATK Haplotype Caller – GATK 4.0.6 https://www.broadinstitute.org/gatk/ 

    Merging GVCFs: Variants are combined using the GATK CombineGVCFs - GATK 4.0.6 and then converted to VCF file by the GATK GenotypeGVCFs - GATK 4.0.6    