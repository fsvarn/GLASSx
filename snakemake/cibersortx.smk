## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CIBERSORTx pipeline for measuring fractions and imputing GEPs
## Author: Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Merge tpms into matrix by batch:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	
rule mergebybatch:
    input:
        expand("results/kallisto/kallisto/aliquot/{aliquot_barcode}/abundance.tsv",aliquot_barcode=manifest.getSelectedAliquots(analyte='R'))
    output:
        protected("results/kallisto/kallisto/final/{batch}_transcript_tpms.tsv")
    log:
        "logs/RNAseq/kallisto/mergebybatch.log"
    message:
        "Merging aliquot TPMs by batch \n"
        "Batch: {wildcards.batch}"
    shell:
    	"""
    	set +o pipefail; 
    	cat {input} | head -1 | sed 's/^/aliquot_barcode\t/' > {output}    	
    	for f in {input}
    	do
    		al=$(echo $f | cut -c 35-64)				#Get aliquot barcode from file name
    		sed "s/^/$al\t/g" $f | tail -n+2 -q >> {output}
    	done 
    	
    	2>{log}
    	"""

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run CIBERSORTx fractions on each batch
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	
rule mergebybatch:
    input:
        "results/kallisto/kallisto/final/{wildcards.batch}_transcript_tpms.tsv"
    output:
        fractions ="results/cibersortx/fractions/CIBERSORTx_{wildcards.batch}_Adjusted.txt",
        mixtures ="results/cibersortx/fractions/CIBERSORTx_{wildcards.batch}_Mixtures_Adjusted.txt"
    params:
    	wd = "/projects/verhaak-lab/GLASS-III/data/subset/cibersortx_fractions_glass/"
    	sigmatrix = "neftel_cibersortx_sig_08032020.txt"
    	bmode = TRUE
    	smode = FALSE	
    	verbose = TRUE
    log:
        "logs/RNAseq/kallisto/cibersortx_fractions_{wildcards.batch}.log"
    message:
        "Running CIBERSORTx fractions module \n"
        "Batch: {wildcards.batch}"
    shell:
    	"""
		module load singularity

		singularity exec \
		-B {params.wd}:/src/data \
		-B /projects/verhaak-lab/GLASS-III/results/cibersortx/fractions/GLASS/:/src/outdir \
		{config.cibersortx_dir}/fractions_latest.sif "/src/CIBERSORTxFractions" \
		--username {config.csx_username} \
		--token  {config.csx_token} \
		--mixture {input} \
		--sigmatrix {params.sigmatrix} \
		--perm 100 \
		--label {wildcards.batch} \
		--rmbatchBmode {params.bmode} \
		--rmbatchSmode {params.smode} \
		--verbose TRUE

    	2>{log}
    	"""

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run CIBERSORTx hi-res mode on each batch
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
	
rule mergebybatch:
    input:
        "results/kallisto/kallisto/final/{wildcards.batch}_transcript_tpms.tsv"
    output:
        fractions ="results/cibersortx/fractions/CIBERSORTx_{wildcards.batch}_Adjusted.txt",
        mixtures ="results/cibersortx/fractions/CIBERSORTx_{wildcards.batch}_Mixtures_Adjusted.txt"
    params:
    	wd = "/projects/verhaak-lab/GLASS-III/data/subset/cibersortx_fractions_glass/"
    	sigmatrix = "neftel_cibersortx_sig_08032020.txt"
    	bmode = TRUE
    	smode = FALSE	
    	verbose = TRUE
    log:
        "logs/RNAseq/kallisto/cibersortx_fractions_{wildcards.batch}.log"
    message:
        "Running CIBERSORTx fractions module \n"
        "Batch: {wildcards.batch}"
    shell:
    	"""
		module load singularity

		singularity exec \
		-B {params.wd}:/src/data \
		-B /projects/verhaak-lab/GLASS-III/results/cibersortx/fractions/GLASS/:/src/outdir \
		{config.cibersortx_dir}/fractions_latest.sif "/src/CIBERSORTxFractions" \
		--username {config.csx_username} \
		--token  {config.csx_token} \
		--mixture {input} \
		--sigmatrix {params.sigmatrix} \
		--perm 100 \
		--label {wildcards.batch} \
		--rmbatchBmode {params.bmode} \
		--rmbatchSmode {params.smode} \
		--verbose TRUE

    	2>{log}
    	"""

