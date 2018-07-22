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

outputs:
    result:
        type: File
        outputSource: picard_addrep/result
steps:
    bwa_mem:
        doc: ALIGNMENT - fq to sam
        run: ../tools/bwa-mem.cwl
        in:
            reference_fa_in: reference_fa_in
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
    
    #picard_fixmateinfo:
    #   doc: SAM/BAM PROCESSING - Fix Mate Information
    #   run: ../tools/picard-fixmateinfo.cwl
    #   in:
    #       sorted_in: picard_sortsam/result
    #       fixmateinfo_out:
    #           valueFrom: $(inputs.sorted_in.nameroot.split("_sorted").slice(0)[0])_fxmt.bam
    #    out: [result]

    picard_markdup:
        doc: SAM/BAM PROCESSING - Mark Duplicates
        run: ../tools/picard-markdup.cwl
        in:
            fixmate_in: picard_sortsam/result
            markdup_out:
                valueFrom: $(inputs.fixmate_in.nameroot.split("_sorted").slice(0)[0])_mkdup.bam
            #fixmate_in: picard_fixmateinfo/result
            #markdup_out:
            #    valueFrom: $(inputs.fixmate_in.nameroot.split("_fxmt").slice(0)[0])_mkdup.bam
        out: [result]

    picard_addrep:
        doc: SAM/BAM PROCESSING - Add Or Replace Read Groups
        run: ../tools/picard-addrep.cwl
        in:
            markdup_in: picard_markdup/result
            addrep_out:
                valueFrom: $(inputs.markdup_in.nameroot.split("_mkdup").slice(0)[0])_addrep.bam
            rgid:
                valueFrom: $(inputs.markdup_in.nameroot.split("_mkdup").slice(0)[0])_id
            rglb:
                valueFrom: $(inputs.markdup_in.nameroot.split("_mkdup").slice(0)[0])_lb
            rgpu:
                valueFrom: $(inputs.markdup_in.nameroot.split("_mkdup").slice(0)[0])_pu
            rgsm:
                valueFrom: $(inputs.markdup_in.nameroot.split("_mkdup").slice(0)[0])_sm
        out: [result]
    
doc: |
    GVCFs Generation is a workflow that creates the gvcf file and to be merge with other gvcfs through Merge GVCFs workflow file.

    1: ALIGNMENT
        a. BWA v0.7.17 --> bwa mem

    2: SAM/BAM PROCESSING
        Picard v2.18.9
        a. Sort 
            --> picard SortSam
        b. Fix Mate Information
            --> picard FixMateInformation
        c. Mark Duplicates
            --> picard MarkDuplicate
        d. Add Or Replace Read Groups
            --> picard AddOrReplaceReadGroups