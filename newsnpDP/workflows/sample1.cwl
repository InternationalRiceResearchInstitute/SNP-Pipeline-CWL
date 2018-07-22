#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
-   class: InlineJavascriptRequirement
-   class: StepInputExpressionRequirement

inputs:
    reference_fa_in: File
    reads1_in: File
    reads2_in: File

outputs:
    result:
        type: File
        outputSource: bwa_mem/result
steps:
    index:
        run: ../indices_generation.cwl
        in:
            reference_fa_in: reference_fa_in
        out: [indices_result]

    gvcfs:
        run: ../tools/bwa-mem.cwl
        in:
            reference_fa_in: reference_fa_in            
            reads1_in: reads1_in
            reads2_in: reads2_in
            sam_out:
                valueFrom: $(inputs.reads1_in.basename.split('.').slice(0)[0])
        out: [result]    

