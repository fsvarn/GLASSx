## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## PRADA pipeline for calling fusions and EGFRvIII variants
## Author: Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run fastp to clean up RNAseq fastq:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule prada_fastp:
    input:
        R1 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[0]),
        R2 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[1])
    output:
        R1 = temp("results/prada/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R1.fastp.fastq"),
        R2 = temp("results/prada/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R2.fastp.fastq"),
        json = "results/prada/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.json",
        html = "results/prada/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.html"
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/prada/fastp/{aliquot_barcode}.{readgroup}.log"
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
## Merge fastq files for input into PRADA
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule prada_merge:
    input:
        R1 = lambda wildcards: expand("results/prada/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R1.fastp.fastq", aliquot_barcode=wildcards.aliquot_barcode, readgroup=manifest.getRGIDs(wildcards.aliquot_barcode)),
        R2 = lambda wildcards: expand("results/prada/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R2.fastp.fastq", aliquot_barcode=wildcards.aliquot_barcode, readgroup=manifest.getRGIDs(wildcards.aliquot_barcode))
    output:
        R1 = temp("results/prada/{aliquot_barcode}/{aliquot_barcode}.end1.fastq"),
    	R2 = temp("results/prada/{aliquot_barcode}/{aliquot_barcode}.end2.fastq")
    log:
        "logs/prada/prada_merge/{aliquot_barcode}.log"
    message:
        "Merging fastq files for input into PRADA \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"""
    	cat {input.R1} > {output.R1} 
    	cat {input.R2} > {output.R2} 
    	2>{log}
    	"""

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run PRADA preprocessing module on fastq files:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule prada_preprocess:
    input:
        R1 = "results/prada/{aliquot_barcode}/{aliquot_barcode}.end1.fastq",
        R2 = "results/prada/{aliquot_barcode}/{aliquot_barcode}.end2.fastq"
    output:
        bam = "results/prada/{aliquot_barcode}/preprocessing/{aliquot_barcode}.withRG.GATKRecalibrated.flagged.bam",
        bai = "results/prada/{aliquot_barcode}/preprocessing/{aliquot_barcode}.withRG.GATKRecalibrated.flagged.bam.bai"
    params:
    	input_dir = "results/prada/{aliquot_barcode}/",
    	sample = "{aliquot_barcode}",
    	output_dir = "results/prada/{aliquot_barcode}/preprocessing/",
        shell = "results/prada/{aliquot_barcode}/preprocessing/script.sh"
    conda:
        "../envs/prada.yaml"
    log:
        "logs/prada/prada_preprocess/{aliquot_barcode}.log"
    message:
        "Running PRADA preprocessing module \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"""
		set +o pipefail; 
    	module load java/1.7.0
    	
    	{config[prada_dir]}/prada-preprocess-bi \
    	-ref {config[prada_conf]} \
		-inputdir {params.input_dir} \
		-sample {params.sample} \
		-tag {params.sample} \
		-step 2_e1_1 \
		-pbs {params.sample} \
		-outdir {params.output_dir} \
		-submit no
		
    	head -1 {params.output_dir}/{params.sample}.pbs > {params.shell}
    	echo 'cd {params.output_dir}' >> {params.shell}
		tail -n +9 {params.output_dir}/{params.sample}.pbs >> {params.shell}
		bash {params.shell}
		
		2>{log}
		"""
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run the fusion module of PRADA to call fusions:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule prada_fusion:
    input:
        "results/prada/{aliquot_barcode}/preprocessing/{aliquot_barcode}.withRG.GATKRecalibrated.flagged.bam"
    output:
        "results/prada/{aliquot_barcode}/fusion/prada.fus.summary.txt"
    params:
    	output_dir = "results/prada/{aliquot_barcode}/fusion/",
    	mismatch = 1,
    	minmapq = 30
    conda:
        "../envs/prada.yaml"
    log:
        "logs/prada/prada_fusion/{aliquot_barcode}.log"
    message:
        "Running PRADA fusion module \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"""
		set +o pipefail; 
    	module load java/1.7.0

    	readlength=$({config[prada_dir]}/tools/samtools-0.1.16/samtools view {input} | \
    	head -n 1000 | \
    	gawk '{{print length($10)}}' | \
    	sort | uniq -c | \
    	perl -ane '$_ =~ s/^[ ]+//g;print $_' | \
    	sort -k 1nr,1nr | \
    	head -1 | \
    	cut -f2 -d " ")
    	
    	junlength=$(python -c "print(int(round($readlength * 0.8)))")

    	{config[prada_dir]}/prada-fusion \
    	-conf {config[prada_conf]} \
    	-bam {input} \
    	-minmapq {params.minmapq} \
    	-mm {params.mismatch} \
    	-junL $junlength \
    	-outdir {params.output_dir}
    	
    	2>{log}
    	"""
    	
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Calculate transcript allele frequency from fusion data:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule prada_taf:
    input:
    	fusion = "results/prada/{aliquot_barcode}/fusion/prada.fus.summary.txt",
        bam = "results/prada/{aliquot_barcode}/preprocessing/{aliquot_barcode}.withRG.GATKRecalibrated.flagged.bam"
    output:
        "results/prada/{aliquot_barcode}/fusion/prada.fus.summary.taf.txt"
    conda:
        "../envs/prada.yaml"
    log:
        "logs/prada/prada_fusion/{aliquot_barcode}.log"
    message:
        "Calculating TAF \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"""
		set +o pipefail; 
    	module load java/1.7.0

		python taf.py \
		--gtf {config[prada_conf]} \
		--fusion {input.fusion}\
		--bam {input.bam} \
		--out {output}
    	
    	2>{log}
    	"""
				
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Run the GUESS-IF module of PRADA to call EGFRvIII variants:
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	
rule prada_guessif:
    input:
        "results/prada/{aliquot_barcode}/preprocessing/{aliquot_barcode}.withRG.GATKRecalibrated.flagged.bam"
    output:
        "results/prada/{aliquot_barcode}/guess-if/{gene}/{gene}.GUESS-IF.summary.txt"
    params:
    	output_dir = "results/prada/{aliquot_barcode}/guess-if/{gene}/",
    	gene = "{gene}",
    	mismatch = 1,
    	minmapq = 30
    conda:
        "../envs/prada.yaml"
    log:
        "logs/prada/prada_guessif/{aliquot_barcode}.{gene}.log"
    message:
        "Running PRADA guess-if module \n"
        "Sample: {wildcards.aliquot_barcode} \n"
        "Gene: {wildcards.gene}"
    shell:
    	"""
		set +o pipefail; 
    	module load java/1.7.0

    	readlength=$({config[prada_dir]}/tools/samtools-0.1.16/samtools view {input} | \
    	head -n 1000 | \
    	gawk '{{print length($10)}}' | \
    	sort | uniq -c | \
    	perl -ane '$_ =~ s/^[ ]+//g;print $_' | \
    	sort -k 1nr,1nr | \
    	head -1 | \
    	cut -f2 -d " ")
    	
    	junlength=$(python -c "print(int(round($readlength * 0.8)))")
    	
    	{config[prada_dir]}/prada-guess-if {params.gene}\
    	-conf {config[prada_conf]} \
    	-inputbam {input} \
    	-minmapq {params.minmapq} \
    	-mm {params.mismatch} \
    	-junL $junlength \
    	-outdir {params.output_dir}
    	
    	2>{log}
    	"""		
