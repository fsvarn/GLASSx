## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## LOHHLA pipeline
## Author: Frederick Varn
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#Defining this variable at the top because LOHHLA uses it in the params and output
MIN_COVERAGE = [5,10,20,30]

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Reformat Optitype output for input into LOHHLA
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule hla_reformat:
    input:
        ancient("results/optitype/HLA_calls/{normal_barcode}/{normal_barcode}_extended_result.tsv")
    output:
        "results/lohhla/prep/{pair_barcode}/{normal_barcode}_lohhla_hlas"
    conda:
        "../envs/lohhla.yaml"
    log:
        "logs/hla_reformat/{pair_barcode}.{normal_barcode}.log"
    message:
        "Reformatting HLA table for LOHHLA input \n"
        "Pair: {wildcards.pair_barcode} \n"
        "Sample: {wildcards.normal_barcode}"
    shell:
        """
        Rscript R/snakemake/hla_optitype2lohhla.R {input} {output}
        2>{log}
        """
    				
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Create solutions.txt file for LOHHLA that contains purity and ploidy (currently using TITAN)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule lohhla_solutions:
    input:
        "results/cnv/titanfinal/params/{pair_barcode}.params.txt",
    output:
        "results/lohhla/prep/{pair_barcode}/{pair_barcode}_solutions.txt"
    params:
    	tumor_barcode = lambda wildcards: "{tumor_barcode}".format(tumor_barcode = manifest.getTumor(wildcards.pair_barcode))
    conda:
        "../envs/lohhla.yaml"
    log:
        "logs/lohhla_solutions/{pair_barcode}.log"
    message:
        "Extracting purity and ploidy info for LOHHLA \n"
        "Sample: {wildcards.pair_barcode}"
    shell:
        """
        Rscript R/snakemake/lohhla_titan_solutions.R {input} {params.tumor_barcode} {output}
        2>{log}
        """					

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Run LOHHLA
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

rule lohhla_run:
    input:
        normal = lambda wildcards: "results/align/bqsr/{normal_barcode}.realn.mdup.bqsr.bam".format(normal_barcode = manifest.getNormal(wildcards.pair_barcode)),
        normal_ind = lambda wildcards: "results/align/bqsr/{normal_barcode}.realn.mdup.bqsr.bai".format(normal_barcode = manifest.getNormal(wildcards.pair_barcode)),
        tumor = lambda wildcards: "results/align/bqsr/{tumor_barcode}.realn.mdup.bqsr.bam".format(tumor_barcode = manifest.getTumor(wildcards.pair_barcode)),
        tumor_ind = lambda wildcards: "results/align/bqsr/{tumor_barcode}.realn.mdup.bqsr.bai".format(tumor_barcode = manifest.getTumor(wildcards.pair_barcode)),
        hlas = lambda wildcards: "results/lohhla/prep/{{pair_barcode}}/{normal_barcode}_lohhla_hlas".format(normal_barcode = manifest.getNormal(wildcards.pair_barcode)),
        solns = "results/lohhla/prep/{pair_barcode}/{pair_barcode}_solutions.txt"
    output:
        "results/lohhla/final/{pair_barcode}/{minCoverageFilter}/{pair_barcode}.{minCoverageFilter}.DNA.HLAlossPrediction_CI.xls"
    params:
    	normal_barcode = lambda wildcards: "{normal_barcode}".format(normal_barcode = manifest.getNormal(wildcards.pair_barcode)),
    	tumor_barcode = lambda wildcards: "{tumor_barcode}".format(tumor_barcode = manifest.getTumor(wildcards.pair_barcode)),
        output_dir = "results/lohhla/final/{pair_barcode}/{minCoverageFilter}/",
        pairID = "{pair_barcode}",
        mappingStep = 'TRUE',
        fishingStep= 'TRUE',
        cleanUp = 'TRUE'
    conda:
        "../envs/lohhla.yaml"
    log:
        "logs/lohhla/{pair_barcode}.{minCoverageFilter}.log"
    message:
        "Running LOHHLA \n"
        "Pair: {wildcards.pair_barcode} \n"
        "Minimum coverage filter: {wildcards.minCoverageFilter}"
    shell:
        """
        mkdir -p {config[lohhla_workdir]}{wildcards.pair_barcode}/{wildcards.minCoverageFilter}
        cp -P {input.normal} {config[lohhla_workdir]}{wildcards.pair_barcode}/{wildcards.minCoverageFilter}/{params.normal_barcode}.bam
        cp -P {input.normal_ind} {config[lohhla_workdir]}{wildcards.pair_barcode}/{wildcards.minCoverageFilter}/{params.normal_barcode}.bai
        cp -P {input.tumor} {config[lohhla_workdir]}{wildcards.pair_barcode}/{wildcards.minCoverageFilter}/{params.tumor_barcode}.bam
        cp -P {input.tumor_ind} {config[lohhla_workdir]}{wildcards.pair_barcode}/{wildcards.minCoverageFilter}/{params.tumor_barcode}.bai
        
        Rscript {config[lohhla_rscript]} \
        --patientId {wildcards.pair_barcode} \
        --outputDir {params.output_dir} \
        --normalBAMfile {config[lohhla_workdir]}{wildcards.pair_barcode}/{wildcards.minCoverageFilter}/{params.normal_barcode}.bam \
        --BAMDir {config[lohhla_workdir]}{wildcards.pair_barcode}/{wildcards.minCoverageFilter}  \
        --hlaPath {input.hlas} \
        --HLAfastaLoc {config[hla_all_fasta]}  \
        --HLAexonLoc {config[hla_exon_loc]} \
        --CopyNumLoc {input.solns} \
        --mappingStep {params.mappingStep} \
        --minCoverageFilter {wildcards.minCoverageFilter} \
        --fishingStep {params.fishingStep} \
        --cleanUp {params.cleanUp} \
        --gatkDir /projects/varnf/SofWar/anaconda3/envs/lohhla/bin/ \
        --novoDir /projects/varnf/SofWar/anaconda3/envs/lohhla/bin/
                        
        2>{log}
        """

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## Format LOHHLA output (doing this externally now)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# 
# rule lohhla_clean:
#     input:
#         expand("results/lohhla/final/{{pair_barcode}}/{minCoverageFilter}/{{pair_barcode}}.{minCoverageFilter}.DNA.HLAlossPrediction_CI.xls", minCoverageFilter=MIN_COVERAGE, allow_missing=True)
#     output:
#         "results/lohhla/final/{pair_barcode}/{pair_barcode}.DNA.HLAlossPrediction_CI.txt"
#     conda:
#         "../envs/lohhla.yaml"
#     log:
#         "logs/lohhla/{pair_barcode}.log"
#     message:
#         "Cleaning LOHHLA output \n"
#         "Pair: {wildcards.pair_barcode} \n"
#     shell:
#         """
#         set +o pipefail; 
#     	cat {input} | head -1 | sed 's/^/pair_barcode\tcoverage_filter\t/' > {output}
#     	for f in {input}
#     	do
#     		cov=$(echo $f | cut -d '/' -f5 | cut -d '.' -f2)
#     		pa="{wildcards.pair_barcode}\t${{cov}}"
#     		sed "s/^/$pa\t/g" $f | tail -n+2 -q >> {output}
#     	done 
#         2>{log}
#         """
