#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
-   class: InlineJavascriptRequirement
-   class: StepInputExpressionRequirement

inputs:
    reference_fa_in: File

outputs:
    indices_result:
        type: File
        outputSource: index_bwa/bwa_indices_result
steps:
    index_fai:
        run: ../tools/samtools-faidx.cwl
        in:
            reference_fa_in: reference_fa_in
        out:
            [fai_result]
    
    index_dict:
        run: ../tools/picard-CSD.cwl
        in:
            reference_fa_in: index_fai/fai_result
            dict_out: 
                valueFrom: $(inputs.reference_fa_in.nameroot).dict 
        out:
            [dict_result]
    
    index_bwa:
        run: ../tools/bwa-index.cwl
        in:
            reference_fa_in: index_dict/dict_result
        out:
            [bwa_indices_result]
    
doc: |
    Indices Generation is a workflow that creates the index files the reference fasta file needed before the invoking the SNP-Discovery-Pipeline.

    1:  Samtools v1.8 --> samtools faidx
    2:  Picard v2.18.9 --> picard CreateSequenceDictionary
    3:  BWA v0.7.17 --> bwa index 