# What is Common Workflow Language or CWL?
Common Workflow Language or CWL is an open standard designed to express workflows and their tooling in groups of YAML structured text files. More details https://www.commonwl.org/

# What is SNP-Pipeline-CWL?
SNP-Pipeline-CWL is the converted form of the IRRI’s SNP-Discovery-Pipeline. It is composed two CWL Pipeline: (1) current SNP pipeline that uses the MERGE BAM step and (2) new SNP pipeline that uses the MERGE GVCF step. Each is used to detect variants including SNPs and indels (insertion/deletions) from next-generation sequencing (NGS) reads. The pipelines includes the following:

1. Indexing: Sam tools is also used to create the some index files – SAM Tools 1.8 https://sourceforge.net/projects/samtools/

2. Alignment: The sequence reads are aligned to the reference genome using BWA (Burrows-Wheeler Aligner) – bwa 0.7.17 https://sourceforge.net/projects/bio-bwa/?source=typ_redirec

3. BAM Processing: A series of intermediate steps to process and prepare the BAM file for variant calling. Picard Tools is used for BAM processing – Picard Tools 2.18.9 https://github.com/broadinstitute/picard 

4. Variant calling: Variants are called using the GATK Unified Genotyper – GATK 4.0.6.0 https://software.broadinstitute.org/gatk/

One of the features of the CWL is interoperability with different environments through the implementation of Docker https://www.docker.com/. Unfortunately this project wasn’t able to implement Docker but the codes can be run in your specific local workstations. 

To run the snp-pipeline-CWL in your personal workstations, you need these things installed in your Ubuntu unit.
1. cwltool - https://github.com/common-workflow-language/cwltool
2. bwa 0.7.17 - https://sourceforge.net/projects/bio-bwa/?source=typ_redirec
3. picard tools 2.18.9 - https://github.com/broadinstitute/picard 
4. samtools 1.8 - https://sourceforge.net/projects/samtools/  
5. gatk 4.0.6.0 - https://software.broadinstitute.org/gatk/

After installing the prerequisites tools and technologies, you may now use these codes in your workstation.

1. Clone repository on your working directory.
	
	`$ git clone https://github.com/jevasQ/internship-project`

2. Go to the directory of the pipeline you want to use.

	a. For using the current SNP Pipeline
	
	`$ cd internship-project/current-SNP-DP`

	b. For using the new SNP Pipeline
	
	`$ cd internship-project/new-SNP-DP`

3. Edit the “input-job” YAML file inside the folder “jobs” and input the directory of your input reference fasta and reads files.

4. Open terminal then invoke the following.

	a. Using the current-SNP-DP for the analysis
	
	`$ cwltool workflows/currentSNPDP.cwl jobs/input-job.yml`

	b. Using the new-SNP-DP for the analysis	
	
	`$ cwltool workflows/newSNPDP.cwl jobs/input-job.yml`
	
5. After the analysis, all the intermediary files created in the process will be deleted and the final output for both of the pipelines is a VCF file.
