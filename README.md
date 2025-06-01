## Author: Sam Blechman
## Date created: January 15th, 2025
## Date last updated: January 17th, 2025

### You just got your sequencing data (.fastq.gz files) back from the Health Sciences Sequencing Core and don't know what to do?? You're in the right place.

This repository houses code for the 7-step sequencing analysis pipeline used by the Wright lab:
1. TRIM: Trim raw sequencing reads based on quality score
2. ASSEMBLE: Make genome assembly for each unique ancestor to be used as reference genome for mutation identification
3. ANNOTATE: Annotate the assembled ancestor genomes
4. BRESEQ: Run breseq to identify mutations that occurred in evolved lineages by mapping their reads to their corresponding annotated and assembled ancestor genomes
5. PARSE_BRESEQ_SUBTRACT_BACKGROUND: Parse breseq output to produce a list of mutations for each evolved lineage + background subtraction (remove mutations that show up in ancestor from evolved lineages)
6. EXTRACT_MUTATION_INFO: Further parse the subtracted breseq output to extract the pertinent mutation information (gene name, position, type of mutation, etc.)
7. IDENTIFY_ENRICHED_MUTATIONS: Identify genes that had mutations at a statistically higher frequency in one group vs. another.

Each step uses different scripts and may be run locally or on the Open Science Grid (OSG or "the grid") for parallelization. Steps that recommend the grid use software that may have environment requirements you won't want to deal with on your own machine.

What is the grid?

The Open Science Grid is a distributed computing network that enables members of our lab to use additional computer resources beyond our personal laptops or servers. Some of the necessary steps in this pipeline require running software using (potentially) hundreds of combinations of inputs (how many evolved lineages do you have?). For example, running breseq on your own laptop for say 200 samples would take several days. However, if each run takes 2 hours, it will be finished in a little over 2 hours on the grid because each of those 200 "jobs" would be completed by a single computer on the grid. This becomes especially helpful if you need to re-run steps after tuning parameters or settings in the software.

To access and use the grid from your personal or work laptop, you must register and get approved. If it's a resource you plan on using beyond this pipeline, consider doing so. However, if this is the only project in the Wright lab for which you require the grid, consider asking a generous lab mate to run your jobs on your behalf. That will require you send them your data and all the necessary scripts personalized to your project and samples.

Each step in more detail, including inputs, outputs, and comments:

STEP 1: TRIM

1. Inputs: raw `[sample_name].fastq.gz` files from HSSC

2. Outputs: 

	a. `[sample_name]_trimmed.fastq.gz` files, name each file according to whether the reads came from an evolved lineage or an ancestor strain (see examples)
	
	b. 'assemble_job_map.txt' containing the file names for the trimmed reads of ancestor strains for assembly (next step)

3. Comments: Run locally, use `R_scripts/001_TrimRawReadsByQuality.R` script


STEP 2: ASSEMBLE

1. Inputs:

	a. trimmed.fastq.gz files (ancestors only)

	b. SPAdes assembly software - `Software/SPAdes.xxx` - there may be a newer version of SPAdes available

	c. Necessary scripts to run this job on the grid - `assemble.sh`, `assemble.sub`, `assemble_job_map.txt` (assemble_job_map.txt as I have it setup right now includes just 1 variable per line: ancestor_trimmed_reads.fastq.gz)

2. Outputs: A directory containing a bunch of stuff. You want to grab the assembled contigs fasta file: `contigs.fasta`. I renamed these files `ancestor_X_contigs.fasta` and put them in the `assemblies/` directory.

3. Comments: Run on grid (could run SPAdes locally, depending on your machine) - there are other assembly programs (just ask Nick) if you feel SPAdes is not ideal or appropriate for your project



STEP 3: ANNOTATE

1. Inputs:

	a. Assembled genome fasta files: `assemblies/ancestor_X_contigs.fasta`

	b. YAML files - not sure how to construct these myself, talk to Nick, I have provided examples and an R script that Nick gave me at some point. They are in the `YAMLfiles/` and `R_scripts/` directories, respectively, but potentially several changes need to be made to this script to update it and make it specific to your organism.

	c. Necessary scripts to run this job on the grid - pgap.sh, pgap.sub, pgap_job_map.txt (`pgap_job_map.txt` as I have it setup right now includes 6 variables per line: job_number, controller.yaml, submol.yaml, output_file_name, input_contigs.fasta)

2. Outputs: Annotated assemblies (`pgap/Annot_X.gbk`)

 	a. I have provided another script called `R_scripts/003_ParseAncestorGBKs.R` which parses the `pgap/Annot_X.gbk` files and grabs a bunch of useful information about the annotated assembled genomes for each strain: number of contigs, contig lengths, gene product names, pgaptmp IDs, direction of transcription, start and end positions, etc. This info is necessary for step 7. The output of `R_scripts/ParseAncestorGBKs.R` is a list of data.frames stored in: `RdataFiles/ancestor_gbk_dataframes.Rdata`.

3. Comments: This is probably the trickiest step that I am least famililar with, ask Nick for help, especially with properly constructing the .yaml files AND ensuring PGAP software and the required environment are good to on the grid.



STEP 4: BRESEQ

1. Inputs: 

	a. `evolved_[sample_name]_trimmed.fastq.gz` files for evolved lineages

	b. `pgap/Annot_X.gbk` files - one annotated and assembled genome per unique ancestor strain, the output from previous steps

	c. Read mapping software - Breseq, Bowtie2 (found in `Software` directory here)

	d. Necessary scripts to run this job on the grid - `map.sh`, `map.sub`, `map_job_map.txt` (`map_job_map.txt` as I have it setup right now incudes 2 variables per line: `evoled_[sample_name]_trimmed.fastq.gz`, `Annot_X.gbk` - I have provided `R_scripts/004_MakeBreseqMapTxt.R` to build the `map_job_map.txt` file programmatically but you will need to modify it to be specific to your project and sample names.)

2. Outputs: Each sample (evolved and ancestor) will have a breseq_output directory that contains a bunch of stuff that breseq/bowtie2 created.

3. Comments: Run on grid



STEP 5: PARSE_BRESEQ_SUBTRACT_BACKGROUND

1. Inputs: breseq_output directory for each sample - I have examples in the `breseq_output` directory, they are pretty big files

2. Outputs: `mutations_subtracted_ancestor/mutations_list.txt` files containing tables of the list of valid mutations for each sample (see comment below)

3. Comments: Run locally, use `R_scripts/005_SubtractAncestorMutations.R` script



STEP 6: EXTRACT_MUTATION_INFO

1. Inputs: `mutations_subtracted_ancestor/[sample_name]_mutations.txt` files

2. Outputs: `RdataFiles/mutations_list.Rdata` file

3. Comments: This script `R_scripts/006_ExtractMutationInfo.R` is a beast. I made it to parse my ~300 samples, which may not have been entirely representative of every combination of outputs that breseq can produce. For help troubleshooting this script if needed, email me at sam.blechman@gmail.com.



STEP 7: IDENTIFY_ENRICHED_MUTATIONS

1. Inputs: `RdataFiles/mutations_list.Rdata` and `RdataFiles/ancestor_gbk_dataframes.Rdata`.

2. Outputs: Two plots showing enriched mutations. You could also save the `data.frame` called `sig_genes` that gets created in this step.

3. Comments: The script `R+scripts/007_IdentifySigGenes.R` looks through the mutations in all the lineages in each group (e.g., VAN-exposed and VAN-unexposed) and counts the number lineages in each group that has any mutation in each gene. It then "scores" each gene according to how enriched each group is in mutations in that gene. To visualize this, two plots are made: i) an X-Y scatter plot showing the proportion of each group that had a mutation in that gene (`mutsXY()`) and ii) a volcano plot showing statistical significance vs. effect size (`volcano()`).


Why is step 5 necessary, you ask? Good question, I answer!

If your experiment is at all similar to mine ([Blechman and Wright, PLoS Pathogens, 2024](https://doi.org/10.1371/journal.ppat.1012422)), then you likely have a handful of ancestor strains, each of which were propagated or evolved in parallel. The ancestor trimmed reads get assembled and annotated, and those assemblies are used to map the corresponding evolved reads against. Interestingly, when you map the ancestor's trimmed reads to the corresponding ancestor assembly (built from those trimmed reads!!), breseq potentially identifies some mutations. These are likely due to errors in library prep, sequencing, or assembly. Any "mutation" that shows up in the ancestor and the evolved lineage, we want to remove it from the evolved lineage. Step 5 does that.


COMMENT ABOUT GRID JOB FILES:
I have provided example executable (.sh), submit (.sub), and mapping (.txt) files for each step that requires jobs to be run on the grid. You need to make your own mapping (.txt) files because your samples and input will be named differently than mine. You also need to check and likely modify the executable (.sh) and submit (.sub) files to change directories, file names, software versions, etc. Please ask for help from someone with experience running this software (Nick, again!).

Good luck!
