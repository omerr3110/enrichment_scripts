#THE FUNCTION RECEIVES A GENE LIST AND COMPARES ITS ENRICHMENT/DEPLETION COMPARED TO THE BACKGROUND (GENOMIC OR USER-SUPPLIED BACKGROUND LIST)

#db: GO_BP, GO_MF, GO_CC, KEGG, OMIM, REACTOME, GAD_Disease
#genelist - the genelist to analyze for enrichment/depletion
#OPTIONAL INPUTS:
#background - a list of genes to compare the genelists to, any gene that is not in the background will be discarded.
#enrORdep - check for enrichment or depletion?
#minGenes - minimum genes found per term. Terms with fewer genes will not be FDR-corrected. If the user chose
#"enr" in enrORdep then the minGenes is the minimum number of *observed* genes from the list found per term.
#If the user chose "dep" in enrORdep then minGenes is the min number of *expected* genes from the list per term.
#outdir - where to save output file. If not entered, default is /Users/dgokhman/HUJI drive/R/R outputs/

TheGrandEnrich_enrich <- function(db,genelist,background,minGenes,outdir,enrORdep)
{
  if (missing(db)) {db <- "GO_BP"}
  if (missing(background)) {background <- c()}
  if (missing(minGenes)) {minGenes <- 0}
  if (missing(outdir)) {outdir <- "/Users/dgokhman/HUJI drive/R/R outputs/"}
  if (missing(enrORdep)) {enrORdep <- "enr"}
  
  #load DB
  geneCol <- 1
  catCol <- 2
  DBdir <- 'Gokhman lab general info/USEFUL DATASETS/Annotation DBs/'
  if (db == "GO_BP") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt",sep="")
  } else if (db == "GO_MF") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2GOTERM_MF_FAT.txt",sep="")
  } else if (db == "GO_CC") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2GOTERM_CC_FAT.txt",sep="")
  } else if (db == "REACTOME-DAVID") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2REACTOME_PATHWAY.txt",sep="")
  } else if (db == "GAD_Disease") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2GAD_DISEASE.txt",sep="")
  } else if (db == "GAD_Disease_class") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2GAD_DISEASE_CLASS.txt",sep="")
  } else if (db == "OMIM") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2OMIM_DISEASE.txt",sep="")
  } else if (db == "KEGG") {
    input_path_gene2category <- paste(DBdir,"KEGG/current/KEGG_pathway_NEW.txt",sep="")
    geneCol <- 2
    catCol <- 5
  }
  if (db %in% c("GO_BP","GO_MF","GO_CC")) {
    gene2cat_World <- read.table(input_path_gene2category, header=F, sep='\t', stringsAsFactors=F, quote="")
  } else {
    gene2cat_World <- read.table(input_path_gene2category, header=T, sep='\t', stringsAsFactors=F, quote="")
  }
  
  #intersect DB with background entered
  if (length(background) > 0) {
    idxBack <- which(toupper(gene2cat_World[,geneCol]) %in% toupper(unlist(background)))
    gene2cat_World <- gene2cat_World[idxBack,]
    #intersectedBackGenes <- intersect(toupper(gene2cat_World[,geneCol]),toupper(unlist(background)))
    #idxBack <- match(toupper(intersectedBackGenes),toupper(gene2cat_World[,geneCol]))
    #gene2cat_World <- gene2cat_World[idxBack,]
  }
  
  #intersect genelist with DB
  uniqueGenes <- unique(unlist(genelist))
  idxGenes <- which(unlist(gene2cat_World[,geneCol]) %in% uniqueGenes)
  gene2cat_intersectedWithGenelist <- gene2cat_World[idxGenes,]
  intersectedGenes <- intersect(unlist(gene2cat_World[,geneCol]),unlist(genelist))
  idxGenes <- match(intersectedGenes,gene2cat_World[,geneCol])
  gene2cat_intersectedWithGenelist_unique <- gene2cat_World[idxGenes,]
  
  colnames <- c("DB","Unique IDs entered","Genes from genelist in DB","Term","Enrichment","Observed","Expected","Genes","P-value","FDR")
  uniqueTerms <- unique(unlist(gene2cat_World[,catCol]))
  uniqueGenes_inDB <- unique(unlist(gene2cat_World[,geneCol]))
  outTable <- as.data.frame(matrix(NA,nrow=length(uniqueTerms),ncol=length(colnames)))
  colnames(outTable) <- colnames
  
  #go over terms
  ListRatio <- nrow(gene2cat_intersectedWithGenelist)/nrow(gene2cat_World)
  for (ll in 1:length(uniqueTerms)) {
    term <- uniqueTerms[ll]
    outTable[ll,"DB"] <- db
    outTable[ll,"Term"] <- term
    outTable[ll,"Unique IDs entered"] <- length(uniqueGenes)
    outTable[ll,"Genes from genelist in DB"] <- nrow(gene2cat_intersectedWithGenelist_unique)
    idxW <- which(gene2cat_World[,catCol] == term) # num of genes linked to term in the DB
    idxO <- which(gene2cat_intersectedWithGenelist[,catCol] == term) # num of genes linked to term in the user-entered list
    outTable[ll,"Observed"] <- length(idxO)
    outTable[ll,"Genes"] <- paste(gene2cat_intersectedWithGenelist[idxO,geneCol],collapse=",") # genes linked to term in the user-entered list
    numGenes4term_DB <- length(idxW)
    if (length(idxW) > 0) {
      outTable[ll,"Expected"] <- ListRatio * numGenes4term_DB
      outTable[ll,"Enrichment"] <- outTable[ll,"Observed"]/outTable[ll,"Expected"]
      if (!is.na(outTable[ll,"Enrichment"])) {
          #print(outTable[ll,"Term"])
          #print(c(outTable[ll,"Observed"],numGenes4term_DB,(length(uniqueGenes_inDB)-numGenes4term_DB),outTable[ll,"Genes from genelist in DB"]))
        
        mat_grand <- matrix(c(outTable[ll,"Observed"], numGenes4term_DB - outTable[ll,"Observed"],
                              outTable[ll,"Genes from genelist in DB"] - outTable[ll,"Observed"],
                            (length(uniqueGenes_inDB)-numGenes4term_DB) - (outTable[ll,"Genes from genelist in DB"] - outTable[ll,"Observed"])),
                          nrow=2, byrow=T)
        if (enrORdep == "enr") {
          outTable[ll,"P-value"] <- stats::fisher.test(mat_grand, alternative = "g")$p.value
        } else if (enrORdep == "dep") {
          outTable[ll,"P-value"] <- stats::fisher.test(mat_grand, alternative = "l")$p.value
        } else if (enrORdep == "both") {
          outTable[ll,"P-value"] <- stats::fisher.test(mat_grand, alternative = "t")$p.value
        }
        
        # if (enrORdep == "enr") {
        #   outTable[ll,"P-value"] <- phyper(outTable[ll,"Observed"],numGenes4term_DB,(length(uniqueGenes_inDB)-numGenes4term_DB),outTable[ll,"Genes from genelist in DB"],lower.tail=F)
        # } else if (enrORdep == "dep") {
        #   outTable[ll,"P-value"] <- phyper(outTable[ll,"Observed"]-1,numGenes4term_DB,(length(uniqueGenes_inDB)-numGenes4term_DB),outTable[ll,"Genes from genelist in DB"],lower.tail=T)
        # } else if (enrORdep == "both"){
        #   if (outTable[ll,"Enrichment"] > 1){
        #     outTable[ll,"P-value"] <- 2*phyper(outTable[ll,"Observed"],numGenes4term_DB,(length(uniqueGenes_inDB)-numGenes4term_DB),outTable[ll,"Genes from genelist in DB"],lower.tail=F)
        #   } else {
        #     outTable[ll,"P-value"] <- 2*phyper(outTable[ll,"Observed"]-1,numGenes4term_DB,(length(uniqueGenes_inDB)-numGenes4term_DB),outTable[ll,"Genes from genelist in DB"],lower.tail=T)
        #   }
        # }
        
      }
      
    }
  }
  
  
  #fdr
  pvals_temp <- outTable[,"P-value"]
  #idxminGenes <- which((outTable[,"Enrichment"] > 1 & outTable[,"Observed"] <= minGenes) | (outTable[,"Enrichment"] < 1 & outTable[,"Expected"] <= minGenes))
  if (enrORdep == "enr") {
    idxminGenes <- which(outTable[,"Observed"] <= minGenes)
  } else if (enrORdep == "dep") {
    idxminGenes <- which(outTable[,"Expected"] <= minGenes)
  } else if (enrORdep == "both"){
    idxminGenes <- which((outTable$Enrichment > 1 & outTable$Observed <= minGenes) | (outTable$Enrichment < 1 & outTable$Expected <= minGenes))
  }
  pvals_temp[idxminGenes] <- NA
  outTable[,"FDR"] <- p.adjust(pvals_temp,method="BH")
  

  
  #WRITE FILE
  if (substr(outdir,nchar(outdir)-1,nchar(outdir)) != "/") {outdir <- paste(outdir,"/",sep="")}
  write.table(outTable, file=paste(outdir,db,"_TheGrandEnrich_enrich.txt",sep=""), 
              sep="\t", quote = F,na = "", row.names = F, col.names = T)
  
  return(outTable)
}
