#Nicole E Soltis

#-------------------------------------------------------
#location of all GEMMA runs for BcAtRNAseq
#then extract relevant effects sizes etc. for trans hotspot analysis!

#At transcripts GEMMA
#/media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachBc_At/output
#kmat script: 
#script: /media/nesoltis/Soltis_AtBc_eQTL/BcAt_RNAGWAS/GEMMA_eachBc_At/A04_runGEMMA_kmat_3genos.sh

rm(list=ls())
setwd("/media/nesoltis/Soltis_npr1coi1_eQTL/data/AtTranscript")

#read in individual GEMMA output files (1 per geno)
#read in all files in folder by pattern matching
#my.files <- list.files(pattern = ".assoc.txt")
#DO FOR COL0 LATER    
  #each phenotype
mytime <- Sys.time()
for (j in c("coi1","npr1")){
  myPhenos <- read.table(paste(j,"/02_GEMMA/binMAF20NA10.fam", sep=""))
  nameslist <- myPhenos[1,6:length(myPhenos)] 
  nameslist[1,1:10] #yes
  myGenes <- as.data.frame(t(nameslist))
  names(myGenes)[1] <- "Gene"
  myGenes$pheno <- 1:nrow(myGenes)
 for (i in 1:23947){
    #actually: 1:23947 for At
    Sys.time()
    my_gemma <- read.table(paste(j,"/04_GEMMAout/",j,"_MAF20NA10_",i,".assoc.txt", sep=""), header=TRUE)
    Sys.time()
    #takes 4 seconds to read 1 phenotype
    #times 24000 = 27 hours
    #take top 1, 10, 100 SNP/phenotype
    my_gemma$pheno <- i
    my_gemma <- merge(my_gemma, myGenes, by = "pheno")
    phenovector <- as.data.frame(my_gemma[,"p_score"])
    names(phenovector)[1] <- paste0(my_gemma[1,"Gene"])
    phenovector <- t(phenovector)
    snpinfo <- t(as.data.frame(my_gemma[,c("chr","ps")]))
    my_gemma_top100 <- my_gemma[order(my_gemma$p_score),]
    my_gemma_top100 <- my_gemma_top100[1:100,]
    row.names(my_gemma_top100) <- c(1:100) + (100*(i-1))
    my_gemma_top10 <- my_gemma_top100[1:10,]
    my_gemma_top1 <- my_gemma_top10[1,]
    #this gives an error but it's fine
    try(ifelse( i == 1, write.table(my_gemma_top100, paste(j,"/05_GEMMAsumm/GEMMA_top100SNPsample.csv", sep=""), sep = ",", col.names = TRUE, row.names=FALSE), write.table(my_gemma_top100, paste(j,"/05_GEMMAsumm/GEMMA_top100SNPsample.csv", sep=""), sep = ",", col.names = FALSE, row.names=FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma_top10, paste(j,"/05_GEMMAsumm/GEMMA_top10SNPsample.csv", sep=""), sep = ",", col.names = TRUE, row.names=FALSE), write.table(my_gemma_top10, paste(j,"/05_GEMMAsumm/GEMMA_top10SNPsample.csv", sep=""), sep = ",", col.names = FALSE, row.names=FALSE, append = TRUE)))
    try(ifelse( i == 1, write.table(my_gemma_top1, paste(j,"/05_GEMMAsumm/GEMMA_top1SNPsample.csv", sep=""), sep = ",", col.names = TRUE, row.names=FALSE), write.table(my_gemma_top1, paste(j,"/05_GEMMAsumm/GEMMA_top1SNPsample.csv", sep=""), sep = ",", col.names = FALSE, row.names=FALSE, append = TRUE)))
    Sys.time()
    try(if( i == 1){
      write.table(snpinfo, paste(j,"/05_GEMMAsumm/GEMMA_allpval_matrix.csv", sep=""), sep=",", col.names=TRUE)
      write.table(phenovector, paste(j,"/05_GEMMAsumm/GEMMA_allpval_matrix.csv", sep=""), sep=",", col.names=FALSE, append=TRUE)
    } else {
      write.table(phenovector, paste(j,"/05_GEMMAsumm/GEMMA_allpval_matrix.csv", sep=""), sep=",", col.names=FALSE, append=TRUE)
    }
    )
 }
}
mytime
Sys.time()