#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
-   class: InlineJavascriptRequirement
-   class: StepInputExpressionRequirement
-   class: SubworkflowFeatureRequirement

inputs:
    reference_fa_in: File
    reads1_in: File
    reads2_in: File

    markdup_in: File
    variant_list_in: File[]
    merged_variant_in: File

outputs:
    # indices: 
    #     type: File
    #     outputSource: indices_generation/indices_result
    
    sam:    
        type: File
        outputSource: bwa_mem/result

    fixmt:
        type: File
        outputSource: picard_fixmateinfo/result

    # test_result:
    #     type: File
    #     outputSource: merging/vcf_result
steps:
    indices_generation:
        doc: Workflow for creating the index files
        run: indices_generation.cwl
        in:
            reference_fa_in: reference_fa_in
        out: [indices_result]

    bwa_mem:
        doc: ALIGNMENT - fq to sam
        run: ../tools/bwa-mem.cwl
        in:
            reference_fa_in: indices_generation/indices_result
            reads1_in: reads1_in
            reads2_in: reads2_in
            sam_out:
                valueFrom: $(inputs.reads1_in.basename.split('.').slice(0)[0]).sam
        out: [result]

    picard_sortsam:
        doc: SAM/BAM PROCESSING - Sort Sam
        run: ../tools/picard-sortsam.cwl
        in: 
            sam_in: bwa_mem/result
            bam_sorted_out:
                valueFrom: $(inputs.sam_in.nameroot)_sorted.bam 
        out: [result]
    
    picard_fixmateinfo:
        doc: SAM/BAM PROCESSING - Fix Mate Information
        run: ../tools/picard-fixmateinfo.cwl
        in:
            sorted_in: picard_sortsam/result
            fixmateinfo_out:
               valueFrom: $(inputs.sorted_in.nameroot.split("_sorted").slice(0)[0])_fxmt.bam
        out: [result]
    # picard_markdup:
    #     doc: SAM/BAM PROCESSING - Mark Duplicates
    #     run: ../tools/picard-markdup.cwl
    #     in:
    #         fixmate_in: picard_sortsam/result
    #         markdup_out:
    #             valueFrom: $(inputs.fixmate_in.nameroot.split("_sorted").slice(0)[0])_mkdup.bam
    #         #fixmate_in: picard_fixmateinfom/result
    #         #markdup_out:
    #             #valueFrom: $(inputs.fixmate_in.nameroot.split("_fxmt").slice(0)[0])_mkdup.bam
    #     out: [result]

    # picard_addrep:
    #     doc: SAM/BAM PROCESSING - Add Or Replace Read Groups
    #     run: ../tools/picard-addrep.cwl
    #     in:
    #         #markdup_in: markdup_in
    #         markdup_in: picard_markdup/result
    #         addrep_out:
    #             valueFrom: $(inputs.markdup_in.nameroot.split("_mkdup").slice(0)[0])_addrep.bam
    #     out: [result]
    
    # gatk_htc_gvcf:
    #     doc: VARIANT CALLING - Haplotype Caller
    #     run: ../tools/gatk-HtC-gvcf.cwl
    #     in:
    #         addrep_in: picard_addrep/result
    #         reference_fa_in: indices_generation/indices_result
    #         gvcf_out:
    #             valueFrom: $(inputs.addrep_in.nameroot.split("_addrep").slice(0)[0]).g.vcf
    #     out: [result]

    # merge:
    #     run: ../tools/gatk-combinegvcfs.cwl
    #     in:
    #         reference_fa_in: indices_generation/indices_result
    #         variant_list_in: variant_list_in
    #         gvcfs_merge_out:
    #             valueFrom: $(inputs.variant_list_in[0].nameroot.split('_').slice(0)[0])_merged.g.vcf
    #     out: [result]

    # merging:
    #     run: gvcfs_merging.cwl
    #     in: 
    #         reference_fa_in: indices_generation/indices_result
    #         variant_list_in: variant_list_in
    #     out: [vcf_result]

doc: |
    Test workflows to check how the created workflows runs in the system and to create the dependency with other workflows.