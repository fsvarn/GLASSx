## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Immune cell inference suite
## Runs different immune cell infiltration inference tools on GLASS expression data
## Tools included:
##	ESTIMATE
##	BASE
## Author: Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run ImmBase to get immune cell levels for CD8+ T cells, B cells, NK cells, macrophages:
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule immbase:
    input:
        "results/kallisto/kallisto/final/transcript_tpms_all_samples.tsv"
    output:
        "results/immune/base/immbase_out.tsv"
    params:
    	sig_matrix = "data/sigs/base4sigs.txt"
    log:
        "logs/immune/immbase.log"
    message:
        "Running ImmBase"
    shell:
    	"(fastp \
		-i {input.R1} \
		-I {input.R2}  \
		-o {output.R1} \
		-O {output.R2}  \
		-j {output.json} \
		-h {output.html}) 2>{log}"

		