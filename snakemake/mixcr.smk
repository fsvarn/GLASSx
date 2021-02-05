## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## MiXCR pipeline for reassembling T cell receptor repertoires
## Author: Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#Define chains for MiXCR exporting. Currently looking at alpha/beta TCRs and all BCRs
xcr_dict = {'tcr': 'TRA,TRB', 'bcr':'IGH,IGL,IGK'}

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run fastp to clean up RNAseq fastq:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule mixcr_fastp:
    input:
        R1 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[0]),
        R2 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[1])
    output:
        R1 = temp("results/mixcr/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R1.fastp.fastq"),
        R2 = temp("results/mixcr/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R2.fastp.fastq"),
        json = "results/mixcr/{aliquot_barcode}/{aliquot_barcode}.{readgroup}_fastp.json",
        html = "results/mixcr/{aliquot_barcode}/{aliquot_barcode}.{readgroup}_fastp.html"
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/RNAseq/fastp/{aliquot_barcode}.{readgroup}.log"
    message:
        "Cleaning up fastq with fastp \n"
        "Sample: {wildcards.aliquot_barcode}.{wildcards.readgroup}"
    shell:
    	"(fastp \
		-i {input.R1} \
		-I {input.R2}  \
		-o {output.R1} \
		-O {output.R2}  \
		-j {output.json} \
		-h {output.html}) 2>{log}"
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge fastq files for input into MiXCR
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule merge_fastq:
    input:
        R1 = lambda wildcards: expand("results/mixcr/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R1.fastp.fastq", aliquot_barcode=wildcards.aliquot_barcode, readgroup=manifest.getRGIDs(wildcards.aliquot_barcode)),
        R2 = lambda wildcards: expand("results/mixcr/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R2.fastp.fastq", aliquot_barcode=wildcards.aliquot_barcode, readgroup=manifest.getRGIDs(wildcards.aliquot_barcode))
    output:
        R1 = temp("results/mixcr/{aliquot_barcode}/{aliquot_barcode}_merge.R1.fastq"),
    	R2 = temp("results/mixcr/{aliquot_barcode}/{aliquot_barcode}_merge.R2.fastq")
    log:
        "logs/mixcr/merge_fastq/{aliquot_barcode}.log"
    message:
        "Merging fastq files for input into MiXCR \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"""
    	cat {input.R1} > {output.R1} 
    	cat {input.R2} > {output.R2} 
    	2>{log}
    	"""

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Align sequencing reads against reference V, D, J, and C genes
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule mixcr_align:
    input:
        R1 = "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_merge.R1.fastq",
        R2 = "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_merge.R2.fastq"
    output:
        "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_alignments.vdjca"
    conda:
        "../envs/mixcr.yaml"
    log:
        "logs/mixcr/align/{aliquot_barcode}.log"
    message:
        "Aligning reads for MiXCR \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"({config[mixcr]} align \
    	-p rna-seq \
    	-s hsa \
    	-OallowPartialAlignments=true \
		{input.R1} \
		{input.R2}  \
		{output}) 2>{log}"
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Assemble partial alignments
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule assemblePartial:
    input:
        "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_alignments.vdjca"
    output:
        r1 = "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_alignments_rescued_1.vdjca",
        r2 = "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_alignments_rescued_2.vdjca"
    conda:
        "../envs/mixcr.yaml"
    log:
        "logs/mixcr/assemblePartial/{aliquot_barcode}.log"
    message:
        "Performing two rounds of partial assembly \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"""
    	{config[mixcr]} assemblePartial \
		{input} \
		{output.r1}
		
		{config[mixcr]} assemblePartial \
		{output.r1} \
		{output.r2}
		2>{log}
		"""
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Extend incomplete TCR CDR3s using germline sequences:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	
rule extend:
    input:
        "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_alignments_rescued_2.vdjca"
    output:
        "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_alignments_rescued_2_extended.vdjca"
    conda:
        "../envs/mixcr.yaml"
    log:
        "logs/mixcr/extend/{aliquot_barcode}.log"
    message:
        "Extending incomplete TCR CDR3s \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"({config[mixcr]} extend \
    	{input} \
    	{output}) 2>{log}"
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Assemble clonotypes:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	
rule assemble:
    input:
        "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_alignments_rescued_2_extended.vdjca"
    output:
        "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_clones.clns"
    conda:
        "../envs/mixcr.yaml"
    log:
        "logs/mixcr/assemble/{aliquot_barcode}.log"
    message:
        "Assembling clonotypes \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"({config[mixcr]} assemble \
    	{input} \
    	{output}) 2>{log}"
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Export TCR and BCR clonotypes separately:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	
rule exportClones:
    input:
        "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_clones.clns"
    output:
        "results/mixcr/{aliquot_barcode}/{aliquot_barcode}_{xcr}_clones.txt"
    params:
    	chains = lambda wildcards: xcr_dict[wildcards.xcr]		#Use dictionary here
    conda:
        "../envs/mixcr.yaml"
    log:
        "logs/mixcr/export/{aliquot_barcode}_{xcr}.log"
    message:
        "Exporting clonotypes \n"
        "Sample: {wildcards.aliquot_barcode}\n"
        "Receptor: {wildcards.xcr}"
    shell:
    	"({config[mixcr]} exportClones \
    	-c {params.chains}\
    	{input} \
    	{output}) 2>{log}"

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge TCR and BCR information and upload to db:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule merge_xcr:
    input:
    	expand("results/mixcr/{aliquot_barcode}/{aliquot_barcode}_{{xcr}}_clones.txt", aliquot_barcode=manifest.getSelectedAliquots(analyte='R'))
    output:
    	protected("results/mixcr/final/merged_clones_{xcr}.txt")
    params:
    	xcr = lambda wildcards: wildcards.xcr
    log:
    	"logs/mixcr/merge_{xcr}.log"
    message:
    	"Merging receptor information into one file and uploading to db \n"
    shell:
    	"""
    	set +o pipefail; 
    	cat {input} | head -1 | sed 's/^/aliquot_barcode\t/' > {output}    	
    	em=$(printf "NA\\t%.0s" {{1..35}})			#Create empty row to add in cases where there are no TCR/BCRs (35 empty fields)
    	for f in {input}
    	do
    		ln=$(wc -l $f | cut -d " " -f 1)
    		al=$(echo $f | cut -c 46-75)				#Get aliquot barcode from file name
    		if [ $ln -eq "1" ]; then echo "$al""\t""$em" >> {output}; fi
    		if [ $ln -gt "1" ]; then sed "s/^/$al\t/g" $f | tail -n+2 -q >> {output}; fi
    	done 
    	
    	Rscript R/snakemake/xcr2db.R {output} {params.xcr}
    	2>{log}
    	"""
