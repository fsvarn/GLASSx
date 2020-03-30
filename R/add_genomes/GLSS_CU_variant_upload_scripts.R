rm(list=ls())

#Use for extracting the files from a specific directory
myDir1 <- "/projects/varnf/GLASS-III/GLASS-III/results/mutect2/geno2db/"
myinf1 <- dir(myDir1)
myinf1 <- myinf1[grep("info.tsv",myinf1)]
myinf1 <- myinf1[grep("GLSS-CU-P",myinf1)]
myinf1 <- myinf1[-grep("GLSS-CU-P103|GLSS-CU-P104",myinf1)] #Remove samples without normals
myoutf1 <- myinf1
myoutf1 <- gsub(".tsv","",myoutf1)
myinf1 <- paste(myDir1,myinf1,sep="")

myTarDir <-"/projects/varnf/TEMP/122019_new_info_GLSS-PD/"
myoutf1 <- paste(myTarDir, myoutf1,".sp",sep="")

for(i in 1:length(myinf1))
{

	script <- paste("Rscript /projects/varnf/GLASS-III/GLASS-III/R/snakemake/info_upload.R",myinf1[i],sep=" ")
	script <- paste("source activate /projects/varnf/SofWar/anaconda3/envs/r_3.6", script, sep="\n")
	conOut = file(myoutf1[i], "w")
	writeLines(script, conOut)
	close(conOut)
}

mysub = paste(myTarDir, "submit.sp", sep="")
conOut = file(mysub, "w")
curLine = rep("sleep 0.2",2*length(myinf1))
curLine[seq(1,length(curLine),by=2)] = paste("qsub -V -lnodes=1:ppn=1,walltime=24:00:00 ", myoutf1, sep="")
writeLines(curLine, conOut)
close(conOut)

comm = paste("chmod u+x ", mysub, sep="")
system(comm)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

rm(list=ls())

#Use for extracting the files from a specific directory
myDir1 <- "/projects/varnf/GLASS-III/GLASS-III/results/mutect2/geno2db/"
myinf1 <- dir(myDir1)
myinf1 <- myinf1[grep("geno.tsv",myinf1)]
myinf1 <- myinf1[grep("GLSS-CU-P",myinf1)]
myinf1 <- myinf1[-grep("GLSS-CU-P103|GLSS-CU-P104",myinf1)] #Remove samples without normals
myoutf1 <- myinf1
myoutf1 <- gsub(".tsv","",myoutf1)
myinf1 <- paste(myDir1,myinf1,sep="")

myTarDir <-"/projects/varnf/TEMP/122019_new_geno_GLSS-PD/"
myoutf1 <- paste(myTarDir, myoutf1,".sp",sep="")

for(i in 1:length(myinf1))
{

	script <- paste("Rscript /projects/varnf/GLASS-III/GLASS-III/R/snakemake/geno_upload.R",myinf1[i],sep=" ")
	script <- paste("source activate /projects/varnf/SofWar/anaconda3/envs/r_3.6", script, sep="\n")
	conOut = file(myoutf1[i], "w")
	writeLines(script, conOut)
	close(conOut)
}

mysub = paste(myTarDir, "submit.sp", sep="")
conOut = file(mysub, "w")
curLine = rep("sleep 0.2",2*length(myinf1))
curLine[seq(1,length(curLine),by=2)] = paste("qsub -V -lnodes=1:ppn=1,walltime=24:00:00 ", myoutf1, sep="")
writeLines(curLine, conOut)
close(conOut)

comm = paste("chmod u+x ", mysub, sep="")
system(comm)
