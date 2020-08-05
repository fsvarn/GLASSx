## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Kallisto pipeline for quantifying noncoding transcript expression in TPM
## Author: Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run fastp to clean up RNAseq fastq:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule fastp_nc:
    input:
        R1 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[0]),
        R2 = lambda wildcards: ancient(manifest.getFASTQ(wildcards.aliquot_barcode, wildcards.readgroup)[1])
    output:
        R1 = temp("results/kallisto/fastp/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R1.fastp.fastq"),
        R2 = temp("results/kallisto/fastp/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.R2.fastp.fastq"),
        json = "results/kallisto/fastp/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.json",
        html = "results/kallisto/fastp/{aliquot_barcode}/{aliquot_barcode}.{readgroup}.html"
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
## Run kallisto to quantify transcripts:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule kallisto_nc:
    input:
        lambda wildcards: expand("results/kallisto/fastp/{aliquot}/{aliquot}.{readgroup}.{read}.fastp.fastq", aliquot=wildcards.aliquot_barcode, readgroup=manifest.getRGIDs(wildcards.aliquot_barcode),read=['R1','R2'])
    output:
        tsv = "results/kallisto/noncoding/aliquot/{aliquot_barcode}/abundance.tsv",
        h5 = "results/kallisto/noncoding/aliquot/{aliquot_barcode}/abundance.h5",
        json = "results/kallisto/noncoding/aliquot/{aliquot_barcode}/run_info.json"
    params:
    	output_dir = "results/kallisto/noncoding/aliquot/{aliquot_barcode}/"
    conda:
        "../envs/kallisto.yaml"
    log:
        "logs/RNAseq/kallisto/aliquot/{aliquot_barcode}.log"
    message:
        "Quantifying transcript with kallisto \n"
        "Sample: {wildcards.aliquot_barcode}"
    shell:
    	"(kallisto quant \
		-i {config[kallisto_nc_idx]} \
		-o {params.output_dir} \
		{input}) 2>{log}"
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge tpms together and upload to db:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	
rule mergetpm_nc:
    input:
        expand("results/kallisto/noncoding/aliquot/{aliquot_barcode}/abundance.tsv",aliquot_barcode=manifest.getSelectedAliquots(analyte='R'))
    output:
        protected("results/kallisto/noncoding/final/nc_transcript_tpms_all_samples.tsv")
    log:
        "logs/RNAseq/kallisto/mergetpm.log"
    message:
        "Merging aliquot TPMs into one file and uploading to database"
    shell:
    	"""
    	set +o pipefail; 
    	cat {input} | head -1 | sed 's/^/aliquot_barcode\t/' > {output}    	
    	for f in {input}
    	do
    		al=$(echo $f | cut -c 36-65)				#Get aliquot barcode from file name
    		sed "s/^/$al\t/g" $f | tail -n+2 -q >> {output}
    	done 
    	
    	Rscript R/snakemake/nctpm2db.R {output}
    	2>{log}
    	"""
		
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Create a tpm matrix to store locally:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	
rule tpm_matrix_nc:
    input:
        expand("results/kallisto/noncoding/aliquot/{aliquot_barcode}/abundance.tsv",aliquot_barcode=manifest.getSelectedAliquots(analyte='R'))
    output:
        transcript = protected("results/kallisto/noncoding/final/nc_transcript_tpm_matrix_all_samples.tsv"),
        gene = protected("results/kallisto/noncoding/final/noncoding_tpm_matrix_all_samples.tsv")
    params:
    	dir_str = "results\/kallisto\/noncoding\/aliquot\/",
    	head_tmp = "results/kallisto/noncoding/final/header.tsv",
    	mat_tmp = "results/kallisto/noncoding/final/transcript_tpm_matrix_all_samples.tsv2"
    log:
        "logs/RNAseq/post/tpm_matrix.log"
    message:
        "Merging aliquot TPMs into a transcript and gene-level matrix"
    shell:
    	"""
		num=$(ls -f1 {input} | wc -l)
		upper=$(echo "$((5 * $num))")
		myseq=$(seq 5 5 $upper | sed 's/^\|$//g' | paste -sd,)
		myseq=$(echo "1,2,"$myseq)
		paste {input} | cut -f $myseq > {output.transcript}

		ls -f1 {input} | sed 's/{params.dir_str}//g' | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){{print "\t$1"}}' | perl -ne 'print "target_id\tlength$_\n"' > {params.head_tmp}
		cat {params.head_tmp} {output.transcript} | grep -v "tpm" > {params.mat_tmp}
		mv {params.mat_tmp} {output.transcript}
		rm -f {params.head_tmp} 
		
		Rscript R/snakemake/transcript2gene.R {output.transcript} {config[ensembl_transcript_mapping]} {output.gene}
		2>{log}
    	"""	
	