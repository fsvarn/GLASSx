rm(list=ls())
mytemplate <- "/projects/verhaak-lab/GLASS-III/bin/seqkit_cleanup.sh"
myDir1 <- "/fastscratch/varnf/tmp_storage/CU_preclean/"
myTarDir <- "/projects/varnf/TEMP/102320_CU_preclean/"
outputDir <- "/fastscratch/varnf/tmp_storage/CU/"

# Original files
mytag <- dir(myDir1)
input <- paste(myDir1, mytag, sep = "")
output <- paste(outputDir, mytag, sep="")
output <- gsub(".gz","", output)
readnum <- gsub(".fq.gz","",sapply(strsplit(mytag,"_"),function(x)x[2]))

conIn = file(mytemplate, "r")
rawcode = readLines(conIn)
close(conIn)

for(i in 1:length(input))
{
	tmp = rawcode
	tmp = gsub("\\{input\\}", input[i], tmp)
	tmp = gsub("\\{output\\}", output[i], tmp)
	tmp = gsub("\\{readnum\\}", readnum[i], tmp)

	data = tmp
	myoutf1 = paste(myTarDir,"job", i, ".sp", sep="")
	conOut = file(myoutf1, "w")
	writeLines(data, conOut)
	close(conOut)
}

mysub = paste(myTarDir, "submit.sp", sep="")
conOut = file(mysub, "w")
curLine = rep("sleep 0.2",2*length(input))
curLine[seq(1,length(curLine),by=2)] = paste("sbatch --nodes=1 --cpus-per-task=1 --time=12:00:00 --mem 12gb job",1:length(input), ".sp", sep="")
writeLines(curLine, conOut)
close(conOut)

comm = paste("chmod u+x ", mysub, sep="")
system(comm)