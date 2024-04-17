#THE FUNCTION RECEIVES A GENE LIST AND COMPARES ITS ENRICHMENT/DEPLETION COMPARED TO THE BACKGROUND (GENOMIC OR USER-SUPPLIED BACKGROUND LIST)

#genelist - the genelist to analyze for enrichment/depletion
#DBversion - 12, 13, 14. Recommended: 14
#curationLevel - confident or confident+tentative. confident: observed in HPO. confident+tentative: observed in HPO or DisGeNET (which includes also mouse and rat phenotypes)
#OPTIONAL INPUTS:
#background - a list of genes to compare the genelists to, any gene that is not in the background will be discarded.
#FreqOfPheno - typical or all. typical: observed in >50% of patients. all: observed at any frequency.
#minGenes - minimum genes per organ. Organs with fewer genes will not be FDR-corrected.
#onlySkel - T or F. T - analyze only skeleton-related organs.
#enrORdep - check for enrichment or depletion?

#set working directory to the Dropbox folder
setwd("C:/Users/nachshoe/Dropbox (Weizmann Institute)")

ORGANizer_Organize <- function(genelist,background,outpath,curationLevel,FreqOfPheno,DBversion,minGenes,onlySkel,enrORdep)
{
  if (missing(background)) {background <- c()}
  if (missing(minGenes)) {minGenes <- 0}
  if (missing(onlySkel)) {onlySkel <- F}
  if (missing(curationLevel)) {curationLevel <- "confident"}
  if (missing(FreqOfPheno)) {FreqOfPheno <- "all"}
  if (missing(DBversion)) {DBversion <- 14}
  if (missing(enrORdep)) {enrORdep <- "both"}
  if(missing(outpath)) {outpath <- "/Users/dgokhman/HUJI drive/R/R outputs/"}
  
  #version
  if (as.numeric(DBversion) == 12) {DBversion <- "12_1"
  } else if (as.numeric(DBversion) == 13) {DBversion <- "13"
  } else if (as.numeric(DBversion) == 14) {DBversion <- "14"
  } else {print("wrong DB version")}
  
  #curationLevel and freq of pheno
  if (tolower(curationLevel) == "confident") {
    if (tolower(FreqOfPheno) == "typical") {
      appendix <- "TYP_v"
    } else if (tolower(FreqOfPheno) == "all") {
      appendix <- "ALL_v"
    } else {
      print("wrong frequency of phenotype entered")
    }
  } else if (tolower(curationLevel) == "confident+tentative") {
    appendix <- "DisHPO_v"
  } else {
    print("wrong curation level entered")
  }
  
  #load ORGANizer DB #these files are just excel versions of the ORGANizer_World mat files
  library("readxl")
  ORGANizer_World_path <- paste("Gokhman lab general info/USEFUL DATASETS/Annotation DBs/Gene ORGANizer/ORGANizer_World_",appendix,DBversion,".xlsx",sep="")
  Organs <- read_excel(ORGANizer_World_path,sheet="Organs",col_names=T)
  Systems <- read_excel(ORGANizer_World_path,sheet="Systems",col_names=T)
  GermLayers <- read_excel(ORGANizer_World_path,sheet="GermLayers",col_names=T)
  Regions <- read_excel(ORGANizer_World_path,sheet="Regions",col_names=T)
  
  #Take only skeletal organs
  if (onlySkel) {
    #this option takes genes that are associated with the skeleton
    #idx <- which(Systems[,"skeleton"] == 1)
    #Organs <- Organs[idx,]
    #Regions <- Regions[idx,]
    #GermLayers <- GermLayers[idx,]
    #Systems <- Systems[idx,]
    #this option takes only skeletal organs
    Organs_skel <- read_excel("Gokhman lab general info/USEFUL DATASETS/Annotation DBs/Gene ORGANizer/ORGANizer_Systems2BodyParts.xlsx",skip=0,col_names=T)
    Organs_skel <- unlist(Organs_skel[which(Organs_skel[,2]==1),1])
    idx <- c(1,2,match(Organs_skel,colnames(Organs)))
    Organs <- Organs[,idx]
  }
  
  
  #intersect background with ORGANizer DB
  if (length(background) > 0) {
    intersectedBackGenes <- intersect(Organs$gene_symbol,unlist(background))
    idxBack <- match(intersectedBackGenes,Organs$gene_symbol)
    Organs <- Organs[idxBack,]
    Systems <- Systems[idxBack,]
    GermLayers <- GermLayers[idxBack,]
    Regions <- Regions[idxBack,]
  }
  
  #intersect genelist with ORGANizer DB
  uniqueGenes <- unique(unlist(genelist))
  intersectedGenes <- intersect(Organs$gene_symbol,unlist(genelist))
  idxGenes <- match(intersectedGenes,Organs$gene_symbol)
  Organs_genelist <- Organs[idxGenes,]
  Systems_genelist <- Systems[idxGenes,]
  GermLayers_genelist <- GermLayers[idxGenes,]
  Regions_genelist <- Regions[idxGenes,]
  
  colnames <- c("Unique IDs entered","Genes from genelist in DB","Type","Body Part","Enrichment","Observed","Expected","Genes","P-value","FDR")
  outTable <- as.data.frame(matrix(NA,nrow=(ncol(Organs)-2+ncol(Systems)-2+ncol(Regions)-2+ncol(GermLayers)-2),ncol=length(colnames)))
  colnames(outTable) <- colnames
  
  outTable[,"Type"] <- c(rep("Organ",times=ncol(Organs)-2),rep("System",times=ncol(Systems)-2),rep("Region",times=ncol(Regions)-2),rep("Germ layer",times=ncol(GermLayers)-2))
  outTable[,"Body Part"] <- c(colnames(Organs)[3:ncol(Organs)],colnames(Systems)[3:ncol(Systems)],colnames(Regions)[3:ncol(Regions)],colnames(GermLayers)[3:ncol(GermLayers)])
  #go over body parts
  ll <- 1
  ll_strt <- ll
  if (length(idxGenes) > 0) {
    for (oo in 3:ncol(Organs)) {
      outTable[ll,"Unique IDs entered"] <- length(uniqueGenes)
      outTable[ll,"Genes from genelist in DB"] <- nrow(Organs_genelist[,oo])
      outTable[ll,"Observed"] <- sum(Organs_genelist[,oo])
      outTable[ll,"Genes"] <- paste(unlist(Organs_genelist[which(Organs_genelist[,oo] == 1),"gene_symbol"]),collapse=",")
      outTable[ll,"Expected"] <- (nrow(Organs_genelist[,oo])/nrow(Organs[,oo])) * sum(Organs[,oo])
      outTable[ll,"Enrichment"] <- outTable[ll,"Observed"]/outTable[ll,"Expected"]
      
      mat_org <- matrix(c(outTable[ll,"Observed"], sum(Organs[,oo]) - outTable[ll,"Observed"],
                          nrow(Organs_genelist[,oo]) - outTable[ll,"Observed"],
                          ((nrow(Organs[,oo])-sum(Organs[,oo]))) - (nrow(Organs_genelist[,oo]) - outTable[ll,"Observed"])),
                        nrow=2, byrow=T)
      if (enrORdep == "enr") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_org, alternative = "g")$p.value
      } else if (enrORdep == "dep") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_org, alternative = "l")$p.value
      } else if (enrORdep == "both") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_org, alternative = "t")$p.value
      }
      # print(colnames(Organs %>% select(oo)))
      # print(mat_org)
      # print(c(outTable[ll,"Observed"], sum(Organs[,oo]), (nrow(Organs[,oo])-sum(Organs[,oo])), nrow(Organs_genelist[,oo])))
      # if (enrORdep == "enr") {
      #   outTable[ll,"P-value"] <- phyper(outTable[ll,"Observed"],sum(Organs[,oo]),(nrow(Organs[,oo])-sum(Organs[,oo])),nrow(Organs_genelist[,oo]),lower.tail=F)
      # } else if (enrORdep == "dep") {
      #   outTable[ll,"P-value"] <- phyper(outTable[ll,"Observed"]-1,sum(Organs[,oo]),(nrow(Organs[,oo])-sum(Organs[,oo])),nrow(Organs_genelist[,oo]),lower.tail=T)
      # } else if (enrORdep == "both"){
      #   if (outTable[ll,"Enrichment"] > 1){
      #     outTable[ll,"P-value"] <- 2*phyper(outTable[ll,"Observed"],sum(Organs[,oo]),(nrow(Organs[,oo])-sum(Organs[,oo])),nrow(Organs_genelist[,oo]),lower.tail=F)
      #   } else{
      #     outTable[ll,"P-value"] <- 2*phyper(outTable[ll,"Observed"]-1,sum(Organs[,oo]),(nrow(Organs[,oo])-sum(Organs[,oo])),nrow(Organs_genelist[,oo]),lower.tail=T)
      #   }
      # }
      ll <- ll+1
    }
    #fdr
    pvals_temp <- outTable[ll_strt:(ll-1),"P-value"]
    idxminGenes <- which((outTable[ll_strt:(ll-1),"Enrichment"] > 1 & outTable[ll_strt:(ll-1),"Observed"] <= minGenes) | (outTable[ll_strt:(ll-1),"Enrichment"] < 1 & outTable[ll_strt:(ll-1),"Expected"] <= minGenes))
    pvals_temp[idxminGenes] <- NA
    outTable[ll_strt:(ll-1),"FDR"] <- p.adjust(pvals_temp,method="BH")
    
    ll_strt <- ll
    for (oo in 3:ncol(Systems)) {
      outTable[ll,"Unique IDs entered"] <- length(uniqueGenes)
      outTable[ll,"Genes from genelist in DB"] <- nrow(Organs_genelist[,oo])
      outTable[ll,"Observed"] <- sum(Systems_genelist[,oo])
      outTable[ll,"Genes"] <- paste(unlist(Systems_genelist[which(Systems_genelist[,oo] == 1),"gene_symbol"]),collapse=",")
      outTable[ll,"Expected"] <- (nrow(Systems_genelist[,oo])/nrow(Systems[,oo])) * sum(Systems[,oo])
      outTable[ll,"Enrichment"] <- outTable[ll,"Observed"]/outTable[ll,"Expected"]
      
      mat_sys <- matrix(c(outTable[ll,"Observed"], sum(Systems[,oo]) - outTable[ll,"Observed"],
                          nrow(Systems_genelist[,oo]) - outTable[ll,"Observed"],
                          ((nrow(Systems[,oo])-sum(Systems[,oo]))) - (nrow(Systems_genelist[,oo]) - outTable[ll,"Observed"])),
                        nrow=2, byrow=T)
      if (enrORdep == "enr") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_sys, alternative = "g")$p.value
      } else if (enrORdep == "dep") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_sys, alternative = "l")$p.value
      } else if (enrORdep == "both") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_sys, alternative = "t")$p.value
      }
      
      # if (enrORdep == "enr") {
      #   outTable[ll,"P-value"] <- phyper(outTable[ll,"Observed"],sum(Systems[,oo]),(nrow(Systems[,oo])-sum(Systems[,oo])),nrow(Systems_genelist[,oo]),lower.tail=F)
      # } else if (enrORdep == "dep") {
      #   outTable[ll,"P-value"] <- phyper(outTable[ll,"Observed"]-1,sum(Systems[,oo]),(nrow(Systems[,oo])-sum(Systems[,oo])),nrow(Systems_genelist[,oo]),lower.tail=T)
      # } else if (enrORdep == "both"){
      #   if (outTable[ll,"Enrichment"] > 1){
      #     outTable[ll,"P-value"] <- 2*phyper(outTable[ll,"Observed"],sum(Systems[,oo]),(nrow(Systems[,oo])-sum(Systems[,oo])),nrow(Systems_genelist[,oo]),lower.tail=F)
      #   } else{
      #     outTable[ll,"P-value"] <- 2*phyper(outTable[ll,"Observed"]-1,sum(Systems[,oo]),(nrow(Systems[,oo])-sum(Systems[,oo])),nrow(Systems_genelist[,oo]),lower.tail=T)
      #   }
      # }
      ll <- ll+1
    }
    #fdr
    pvals_temp <- outTable[ll_strt:(ll-1),"P-value"]
    idxminGenes <- which((outTable[ll_strt:(ll-1),"Enrichment"] > 1 & outTable[ll_strt:(ll-1),"Observed"] <= minGenes) | (outTable[ll_strt:(ll-1),"Enrichment"] < 1 & outTable[ll_strt:(ll-1),"Expected"] <= minGenes))
    pvals_temp[idxminGenes] <- NA
    outTable[ll_strt:(ll-1),"FDR"] <- p.adjust(pvals_temp,method="BH")
    
    ll_strt <- ll
    for (oo in 3:ncol(Regions)) {
      outTable[ll,"Unique IDs entered"] <- length(uniqueGenes)
      outTable[ll,"Genes from genelist in DB"] <- nrow(Organs_genelist[,oo])
      outTable[ll,"Observed"] <- sum(Regions_genelist[,oo])
      outTable[ll,"Genes"] <- paste(unlist(Regions_genelist[which(Regions_genelist[,oo] == 1),"gene_symbol"]),collapse=",")
      outTable[ll,"Expected"] <- (nrow(Regions_genelist[,oo])/nrow(Regions[,oo])) * sum(Regions[,oo])
      outTable[ll,"Enrichment"] <- outTable[ll,"Observed"]/outTable[ll,"Expected"]
      
      mat_reg <- matrix(c(outTable[ll,"Observed"], sum(Regions[,oo]) - outTable[ll,"Observed"],
               nrow(Regions_genelist[,oo]) - outTable[ll,"Observed"],
               ((nrow(Regions[,oo])-sum(Regions[,oo]))) - (nrow(Regions_genelist[,oo]) - outTable[ll,"Observed"])),
             nrow=2, byrow=T)
      if (enrORdep == "enr") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_reg, alternative = "g")$p.value
      } else if (enrORdep == "dep") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_reg, alternative = "l")$p.value
      } else if (enrORdep == "both") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_reg, alternative = "t")$p.value
      }
      
      # if (enrORdep == "enr") {
      #   outTable[ll,"P-value"] <- phyper(outTable[ll,"Observed"],sum(Regions[,oo]),(nrow(Regions[,oo])-sum(Regions[,oo])),nrow(Regions_genelist[,oo]),lower.tail=F)
      # } else if (enrORdep == "dep") {
      #   outTable[ll,"P-value"] <- phyper(outTable[ll,"Observed"]-1,sum(Regions[,oo]),(nrow(Regions[,oo])-sum(Regions[,oo])),nrow(Regions_genelist[,oo]),lower.tail=T)
      # } else if (enrORdep == "both"){
      #   if (outTable[ll,"Enrichment"] > 1){
      #     outTable[ll,"P-value"] <- 2*phyper(outTable[ll,"Observed"],sum(Regions[,oo]),(nrow(Regions[,oo])-sum(Regions[,oo])),nrow(Regions_genelist[,oo]),lower.tail=F)
      #   } else{
      #     outTable[ll,"P-value"] <- 2*phyper(outTable[ll,"Observed"]-1,sum(Regions[,oo]),(nrow(Regions[,oo])-sum(Regions[,oo])),nrow(Regions_genelist[,oo]),lower.tail=T)
      #   }
      # }
      ll <- ll+1
    }
    #fdr
    pvals_temp <- outTable[ll_strt:(ll-1),"P-value"]
    idxminGenes <- which((outTable[ll_strt:(ll-1),"Enrichment"] > 1 & outTable[ll_strt:(ll-1),"Observed"] <= minGenes) | (outTable[ll_strt:(ll-1),"Enrichment"] < 1 & outTable[ll_strt:(ll-1),"Expected"] <= minGenes))
    pvals_temp[idxminGenes] <- NA
    outTable[ll_strt:(ll-1),"FDR"] <- p.adjust(pvals_temp,method="BH")
    
    ll_strt <- ll
    for (oo in 3:ncol(GermLayers)) {
      outTable[ll,"Unique IDs entered"] <- length(uniqueGenes)
      outTable[ll,"Genes from genelist in DB"] <- nrow(Organs_genelist[,oo])
      outTable[ll,"Observed"] <- sum(GermLayers_genelist[,oo])
      outTable[ll,"Genes"] <- paste(unlist(GermLayers_genelist[which(GermLayers_genelist[,oo] == 1),"gene_symbol"]),collapse=",")
      outTable[ll,"Expected"] <- (nrow(GermLayers_genelist[,oo])/nrow(GermLayers[,oo])) * sum(GermLayers[,oo])
      outTable[ll,"Enrichment"] <- outTable[ll,"Observed"]/outTable[ll,"Expected"]
      
      mat_gel <- matrix(c(outTable[ll,"Observed"], sum(GermLayers[,oo]) - outTable[ll,"Observed"],
                          nrow(GermLayers_genelist[,oo]) - outTable[ll,"Observed"],
                          ((nrow(GermLayers[,oo])-sum(GermLayers[,oo]))) - (nrow(GermLayers_genelist[,oo]) - outTable[ll,"Observed"])),
                        nrow=2, byrow=T)
      if (enrORdep == "enr") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_gel, alternative = "g")$p.value
      } else if (enrORdep == "dep") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_gel, alternative = "l")$p.value
      } else if (enrORdep == "both") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_gel, alternative = "t")$p.value
      }
      
      # if (enrORdep == "enr") {
      #   outTable[ll,"P-value"] <- phyper(outTable[ll,"Observed"],sum(GermLayers[,oo]),(nrow(GermLayers[,oo])-sum(GermLayers[,oo])),nrow(GermLayers_genelist[,oo]),lower.tail=F)
      # } else if (enrORdep == "dep") {
      #   outTable[ll,"P-value"] <- phyper(outTable[ll,"Observed"]-1,sum(GermLayers[,oo]),(nrow(GermLayers[,oo])-sum(GermLayers[,oo])),nrow(GermLayers_genelist[,oo]),lower.tail=T)
      # } else if (enrORdep == "both") {
      #   if (outTable[ll,"Enrichment"] > 1){
      #     outTable[ll,"P-value"] <- 2*phyper(outTable[ll,"Observed"],sum(GermLayers[,oo]),(nrow(GermLayers[,oo])-sum(GermLayers[,oo])),nrow(GermLayers_genelist[,oo]),lower.tail=F)
      #   } else{
      #     outTable[ll,"P-value"] <- 2*phyper(outTable[ll,"Observed"]-1,sum(GermLayers[,oo]),(nrow(GermLayers[,oo])-sum(GermLayers[,oo])),nrow(GermLayers_genelist[,oo]),lower.tail=T)
      #   }
      # }
      ll <- ll+1
    }
    #fdr
    pvals_temp <- outTable[ll_strt:(ll-1),"P-value"]
    idxminGenes <- which((outTable[ll_strt:(ll-1),"Enrichment"] > 1 & outTable[ll_strt:(ll-1),"Observed"] <= minGenes) | (outTable[ll_strt:(ll-1),"Enrichment"] < 1 & outTable[ll_strt:(ll-1),"Expected"] <= minGenes))
    pvals_temp[idxminGenes] <- NA
    outTable[ll_strt:(ll-1),"FDR"] <- p.adjust(pvals_temp,method="BH")
  }
  
  #WRITE FILE
  write.table(outTable, file=paste0(outpath,"GeneORGANizer_ORGANize.txt"),
              sep="\t", quote = F,na = "", row.names = F, col.names = T)
  
  return(outTable)
}
