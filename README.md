## Author: Sam Blechman
## Date created: January 15th, 2025
## Date last updated: January 17th, 2025

### You just got your sequencing data (.fastq.gz files) back from the Health Sciences Sequencing Core and don't know what to do?? You're in the right place.

This repository houses code for the 6-step sequencing analysis pipeline used by the Wright lab:
1. TRIM: Trim raw sequencing reads based on quality score
2. ASSEMBLE: Make genome assembly for each unique ancestor to be used as reference genome for mutation identification
3. ANNOTATE: Annotate the assembled ancesotr genomes
4. BRESEQ: Run breseq to identify mutations that occurred in evolved lineages by mapping their reads to their corresponding annotated and assembled ancestor genomes
5. PARSE_BRESEQ_SUBTRACT_BACKGROUND: Parse breseq output to produce a list of mutations for each evolved lineage + background subtraction (remove mutations that show up in ancestor from evolved lineages)
6. EXTRACT_MUTATION_INFO: Further parse the subtracted breseq output to extract the pertinent mutation information (gene name, position, type of mutation, etc.)

Each step uses different scripts and may be run locally or on the Open Science Grid (OSG or "the grid") for parallelization. Steps that recommend the grid use software that may have environment requirements you won't want to deal with on your own machine.

What is the grid?

The Open Science Grid is a distributed computing network that enables members of our lab to use additional computer resources beyond our personal laptops or servers. Some of the necessary steps in this pipeline require running software using (potentially) hundreds of combinations of inputs (how many evolved lineages do you have?). For example, running breseq on your own laptop for say 200 samples would take several days. However, if each run takes 2 hours, it will be finished in a little over 2 hours on the grid because each of those 200 "jobs" would be completed by a single computer on the grid. This becomes especially helpful if you need to re-run steps after tuning parameters or settings in the software.

To access and use the grid from your personal or work laptop, you must register and get approved. If it's a resource you plan on using beyond this pipeline, consider doing so. However, if this is the only project in the Wright lab for which you require the grid, consider asking a generous lab mate to run your jobs on your behalf. That will require you send them your data and all the necessary scripts personalized to your project and samples.

Each step in more detail, including inputs, outputs, and comments:
1. TRIM:

	a. Inputs: raw `[sample_name].fastq.gz` files from HSSC

	b. Outputs: 

		i. `[sample_name]_trimmed.fastq.gz` files, name each file according to whether the reads came from an evolved lineage or an ancestor strain (see examples)
		
		ii. 'assemble_job_map.txt' containing the file names for the trimmed reads of ancestor strains for assembly (next step)
	
	c. Comments: Run locally, use `R_scripts/TrimRawReadsByQuality.R` script

2. ASSEMBLE:
	
	a. Inputs:
	
		i. trimmed.fastq.gz files (ancestors only)
	
		ii. SPAdes assembly software - `Software/SPAdes.xxx` - there may be a newer version of SPAdes available
	
		iii. Necessary scripts to run this job on the grid - `assemble.sh`, `assemble.sub`, `assemble_job_map.txt`
	
			1. assemble_job_map.txt as I have it setup right now includes just 1 variable per line: ancestor_trimmed_reads.fastq.gz
	
	b. Outputs: A directory containing a bunch of stuff. You want to grab the assembled contigs fasta file: `contigs.fasta`. I renamed these files `ancestor_X_contigs.fasta` and put them in the `assemblies/` directory.
	
	c. Comments: Run on grid (could run SPAdes locally, depending on your machine) - there are other assembly programs (just ask Nick) if you feel SPAdes is not ideal or appropriate for your project

3. ANNOTATE:
	
	a. Inputs:
	
		i. Assembled genome fasta files: `assemblies/ancestor_X_contigs.fasta`
	
		ii. YAML files - not sure how to construct these myself, talk to Nick, I have provided examples and an R script that Nick gave me at some point. They are in the `YAMLfiles/` and `R_scripts/` directories, respectively, but potentially several changes need to be made to this script to update it and make it specific to your organism.
	
		iii. Necessary scripts to run this job on the grid - pgap.sh, pgap.sub, pgap_job_map.txt
	
			1. `pgap_job_map.txt` as I have it setup right now includes 6 variables per line: job_number, controller.yaml, submol.yaml, output_file_name, input_contigs.fasta
	
	b. Outputs: Annotated assemblies (`pgap/Annot_X.gbk`)
	
	c. Comments: This is probably the trickiest step that I am least famililar with, ask Nick for help, especially with properly constructing the .yaml files AND ensuring PGAP software and the required environment are good to on the grid.

4. BRESEQ:
	
	a. Inputs: 
	
		i. `evolved_[sample_name]_trimmed.fastq.gz` files for evolved lineages
	
		ii. `pgap/Annot_X.gbk` files - one annotated and assembled genome per unique ancestor strain, the output from previous steps
	
		iii. Read mapping software - Breseq, Bowtie2 (found in `Software` directory here)
	
		iv. Necessary scripts to run this job on the grid - `map.sh`, `map.sub`, `map_job_map.txt`
	
			1. `map_job_map.txt` as I have it setup right now incudes 2 variables per line: `evoled_[sample_name]_trimmed.fastq.gz`, `Annot_X.gbk` - I have provided `R_scripts/MakeBreseqMapTxt.R` to build the `map_job_map.txt` file programmatically but you will need to modify it to be specific to your project and sample names.
	
	b. Outputs: Each sample (evolved and ancestor) will have a breseq_output directory that contains a bunch of stuff that breseq/bowtie2 created.
	
	c. Comments: Run on grid

5. PARSE_BRESEQ_SUBTRACT_BACKGROUND:
	
	a. Inputs: breseq_output directory for each sample - I have examples in the `breseq_output` directory, they are pretty big files
	
	b. Outputs: `mutations_subtracted_ancestor/mutations_list.txt` files containing tables of the list of valid mutations for each sample (see comment below)
	
	c. Comments: Run locally, use `R_scripts/SubtractAncestorMutations.R` script

6. EXTRACT_MUTATION_INFO:
	
	a. Inputs: `mutations_subtracted_ancestor/[sample_name]_mutations.txt` files
	
	b. Outputs: `RdataFiles/mutations_list.Rdata` file
	
	c. Comments: This script `R_scripts/ExtractMutationInfo.R` is a beast. I made it to parse my ~300 samples, which may not have been entirely representative of every combination of outputs that breseq can produce. For help troubleshooting this script if needed, email me at sam.blechman@gmail.com.


I have also provided another script called `R_scripts/ParseAncestorGBKs.R` which parses the `pgap/Annot_X.gbk` files and grabs a bunch of useful information for each strain: number of contigs, contig lengths, gene product names, pgaptmp IDs, direction of transcription, start and end positions, etc. This info may be useful when working through the `RdataFiles/mutations_list.Rdata` data.


Why is step 5 necessary, you ask? Good question, I answer!

If your experiment is at all similar to mine ([Blechman and Wright, PLoS Pathogens, 2024](https://doi.org/10.1371/journal.ppat.1012422)), then you likely have a handful of ancestor strains, each of which were propagated or evolved in parallel. The ancestor trimmed reads get assembled and annotated, and those assemblies are used to map the corresponding evolved reads against. Interestingly, when you map the ancestor's trimmed reads to the corresponding ancestor assembly (built from those trimmed reads!!), breseq potentially identifies some mutations. These are likely due to errors in library prep, sequencing, or assembly. Any "mutation" that shows up in the ancestor and the evolved lineage, we want to remove it from the evolved lineage. Step 5 does that.


COMMENT ABOUT GRID JOB FILES:
I have provided example executable (.sh), submit (.sub), and mapping (.txt) files for each step that requires jobs to be run on the grid. You need to make your own mapping (.txt) files because your samples and input will be named differently than mine. You also need to check and likely modify the executable (.sh) and submit (.sub) files to change directories, file names, software versions, etc. Please ask for help from someone with experience running this software (Nick, again!).

Good luck!