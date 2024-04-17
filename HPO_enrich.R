# This function receives a gene list and returns the enrichment and p-value of each HPO associated with the list

# genelist - genes to analyze (enter gene symbols)
#OPTIONAL INPUTS:
# outpuath - dir where the output will be saved
# appendix - any additional string to add to the output file
# background - the background list to compare genelist to. Genes that are not in the background list will be discarded. Default: all genes.
# minGenes - HPO terms linked to fewer genes than this number will not be FDR-corrected.
# TYPofALL - typical or all. typical: observed in >50% of patients. all: observed at any frequency.
# HPO_dbVer - 12, 13 or 14. Recommended: 14.
# bodyparts - a list of bodyparts. When the user supplies a list of bodyparts, the function will only analyze phenotypes associated with these body parts. Any body part that appears in Gene ORGANizer can be used.
# systems - a list of systems. When the user supplies a list of systems, the function will only analyze phenotypes associated with these systems. Any system that appears in Gene ORGANizer can be used.
# directional - T or F. T: analyze only directional phenotypes (i.e., phenotypes that can be described on a scale)
#enrORdep - check for enrichment or depletion?
#locORser - running locally or on the server - affects on DBdir
#comp - string variable - "diff_activity" - diff_active vs. active or "activity" - active vs. all

HPO_enrich <- function(genelist, outpath, appendix, background, minGenes, TYPorALL, HPO_dbVer, systems, bodyparts, directional, enrORdep,locORser,comp)
{
  if (missing(background)) {background <- c()}
  if (missing(appendix)) {appendix <- ""}
  if (missing(minGenes)) {minGenes <- 5}
  if (missing(bodyparts)) {bodyparts <- "all"} else {bodyparts <- tolower(bodyparts)}
  if (missing(systems)) {systems = "all"}
  if (missing(TYPorALL)) {TYPorALL <- "ALL"} else {TYPorALL <- toupper(TYPorALL)}
  if (missing(HPO_dbVer)) {HPO_dbVer <- 14}
  if (missing(directional)) {directional <- F}
  if (missing(outpath)) {outpath <- "/Users/dgokhman/HUJI drive/R/R outputs/"}
  if (missing(enrORdep)) {enrORdep <- "enr"}
  if (missing(locORser)) {locORser <- "loc"}
  if (missing(comp)) {comp <- "diff_activity"}
  genelist <- unique(unlist(genelist))
  background <- unique(unlist(background))
  if (locORser=="ser"){
    load(paste("/home/labs/davidgo/Collaboration/USEFUL DATASETS/Annotation DBs/HPO/HPO_world/HPO_World_",TYPorALL,"_",HPO_dbVer,".RData",sep=""))
  } else{  load(paste("Gokhman lab general info/USEFUL DATASETS/Annotation DBs/HPO/HPO_world/HPO_World_",TYPorALL,"_",HPO_dbVer,".RData",sep=""))
    }
  
  #CONCATENATE BODY PARTS
  bodyparts_conc <- paste0(bodyparts,collapse = ", ")
  systems_conc <- paste0(systems,collapse = ", ")
  
  #REMOVE UNINFORMATIVE HPO TERMS
  hpos2remove <- c("Autosomal recessive inheritance","Autosomal dominant inheritance","X-linked dominant inheritance","X-linked recessive inheritance","Somatic mutation","X-linked inheritance","Polygenic inheritance","Multifactorial inheritance")
  idx2remove <- c()
  for (hpo in hpos2remove) {
    idx <- which(HPO_Phenos == hpo)
    if (length(idx) > 0) {idx2remove <- c(idx2remove,idx)}
  }
  HPO_World <- HPO_World[-idx2remove,]
  HPO_IDs <- HPO_IDs[-idx2remove]
  HPO_Phenos <- HPO_Phenos[-idx2remove]
  
  #LOOK FOR HPOs LINKED TO BODY PARTS
  library("readxl")
  if (HPO_dbVer == 14) {sheet <- "Pheno - build1268 - v14"
  } else if (HPO_dbVer == 13) {sheet <- "Pheno - build115 - v13"
  } else if (HPO_dbVer == 12) {sheet <- "Pheno - v12"
  }
  if (locORser=="ser"){
    HPO_ORGANizer_anno <- read_excel("/home/labs/davidgo/Collaboration/USEFUL DATASETS/Annotation DBs/HPO/specific phenotype ontology.xlsx",skip=0,sheet=sheet,col_names=T)
  } else{  HPO_ORGANizer_anno <- read_excel("Gokhman lab general info/USEFUL DATASETS/Annotation DBs/HPO/specific phenotype ontology.xlsx",skip=0,sheet=sheet,col_names=T)
  }
  HPO_ORGANizer_anno <- HPO_ORGANizer_anno[,c("HPO-ID","system","organ / body part")]
  
  #FILTER TO TAKE ONLY THE USER-ENTERED BODY PART
  hpos2keep_sys <- c()
  if ("all" %in% systems) {
    hpos2keep_sys <- HPO_IDs; print("not filtering based on system")
  } else {
    print(paste0("taking only phenotypes associated with the systems: ",systems))
    for (system in systems) {
      for (row in 1:nrow(HPO_ORGANizer_anno)) {
        pos <- gregexpr(system,HPO_ORGANizer_anno[row,"system"],fixed=T)[[1]]
        if (!is.na(pos) & pos[1] != -1) {
          hpos2keep_sys <- c(hpos2keep_sys,HPO_ORGANizer_anno[row,"HPO-ID"])
        }
      }
    }
  }
  
  hpos2keep <- c()
  if ("all" %in% bodyparts) {
    print("not filtering based on body parts")
    hpos2keep <- HPO_IDs
  }  else {
    print(paste0("taking phenotypes associated with these bodyparts: ",bodyparts_conc))
    print(paste0("taking phenotypes associated with these systems: ",systems_conc))
    for (row in 1:nrow(HPO_ORGANizer_anno)) {
      for (bodypart in bodyparts) {
        pos <- regexpr(bodypart, HPO_ORGANizer_anno[row,"organ / body part"],fixed=T)[[1]]
        if (!is.na(pos) & pos[1] != -1) {
          hpos2keep <- c(hpos2keep,HPO_ORGANizer_anno[row,"HPO-ID"])
        }
      }
    }
  }
  
  hpo2keep_directional <- HPO_IDs
  if (directional) {
    load(paste("DerivedTraits_",TYPorALL,"_",HPO_dbVer,".RData",sep=""))
    hpo2keep_directional <- as.character(unique(gene2hpo2trait[,"HPO ID"]))
  }
  
  hpos2keep <- intersect(intersect(hpos2keep,hpos2keep_sys),hpo2keep_directional)
  bodypartHPOs <- intersect(HPO_IDs,hpos2keep)
  idxhpos2keep <- match(bodypartHPOs,HPO_IDs)
  HPO_IDs <- HPO_IDs[idxhpos2keep]
  HPO_Phenos <- HPO_Phenos[idxhpos2keep]
  HPO_World <- HPO_World[idxhpos2keep,]
  HPO_Genes <- HPO_Genes[idxhpos2keep]
  HPO_ENTREZ_IDs <- HPO_ENTREZ_IDs[idxhpos2keep]
  hposPerGene <- colSums(HPO_World)
  idx <- which(colSums(HPO_World) == 0)
  if (length(idx) > 0) {
    HPO_World <- HPO_World[,-idx]
    HPO_Genes <- HPO_Genes[-idx]
    HPO_ENTREZ_IDs <- HPO_ENTREZ_IDs[-idx]
  }
  
  # CREATE BACKGROUND
  if (length(background) > 0) {
    idx <- c()
    for (gg in 1:length(background)) {
      gene <- toupper(background[gg])
      ii <- which(HPO_Genes == gene)
      if (length(ii) > 0) {
        idx <- rbind(idx,ii)
      }
    }
    idx <- unlist(idx)
    HPO_World <- HPO_World[,idx]
    HPO_Genes <- HPO_Genes[idx]
    HPO_ENTREZ_IDs <- HPO_ENTREZ_IDs[idx]
  }
  
  # CREATE USER-ENTERED SUB-TABLE OF HPO_WORLD WITH GENE-HPO ASSOCIATIONS BASED ON ENTERED GENE LIST
  idx <- c()
  for (gg in 1:length(genelist)) {
    gene <- toupper(genelist[gg])
    ii <- which(HPO_Genes == gene)
    if (length(ii) > 0) {
      idx <- rbind(idx,ii)
    }
  }
  idx <- unlist(idx)
  HPO_Entered <- HPO_World[,idx]
  HPO_Entered_Genes <- HPO_Genes[idx]
  HPO_Entered_ENTREZ_IDs <- HPO_ENTREZ_IDs[idx]
  
  #GO OVER EACH HPO AND COMPUTE P-VALS
  # N_world - number of genes found both within background gene list (active genes) and HPO's gene list
  # N_entered - number of genes found both within test gene list (differentially active) and HPO's gene list
  # n_world - number of genes from background gene list found to be associated with a certain phenotype
  # Observed (n_entered) - number of genes from test gene list found to be associated with a certain phenotype
  
  cols <- c("Pheno","HPO_ID","Enrichment","N_world","n_world","N_entered","Observed (n_entered)","Expected","Genes","P-value","FDR")
  StatMat <- matrix(NA,nrow=nrow(HPO_World),ncol=length(cols))
  colnames(StatMat) <- cols
  StatMat <- as.data.frame(StatMat)
  
  N_world <- ncol(HPO_World)
  N_entered <- ncol(HPO_Entered)
  if (length(N_entered) > 0) {
    for (hpo in 1:nrow(HPO_World)) {
      StatMat$N_world[hpo] <- N_world      
      StatMat$N_entered[hpo] <- N_entered
      StatMat$Pheno[hpo] <- HPO_Phenos[hpo]
      StatMat$HPO_ID[hpo] <- HPO_IDs[hpo]
      n_world <- sum(HPO_World[hpo,])
      StatMat$n_world[hpo] <- n_world
      n_entered <- sum(HPO_Entered[hpo,])
      StatMat$Enrichment[hpo] <- (n_entered/N_entered)/(n_world/N_world)
      StatMat[hpo,"Observed (n_entered)"] <- n_entered
      idx <- which(HPO_Entered[hpo,] == 1)
      StatMat[hpo,"Genes"] <- paste(HPO_Entered_Genes[idx],collapse=",")
      StatMat$Expected[hpo] <- (n_world/N_world)*N_entered  #number of active genes associated to phenotype times all diff active genes over all active genes
      
      mat_hpo <- matrix(c(n_entered, n_world - n_entered, #number of genes associated to the phenotype that are active but not diff active
                          N_entered - n_entered,          #number of diff active genes that aren't associated to the phenotype
                          (N_world-n_world) - (N_entered - n_entered)),  #number of active, non-diff active genes that aren't associated to the phenotype
                        nrow=2, byrow=T)
      row_sums=rowSums(mat_hpo)
      col_sums=colSums(mat_hpo)
      tot_sum=sum(mat_hpo)
      exp_mat=matrix(c(row_sums[1]*col_sums[1]/tot_sum,row_sums[2]*col_sums[1]/tot_sum,
                       row_sums[1]*col_sums[2]/tot_sum,row_sums[2]*col_sums[2]/tot_sum),nrow=2)
      bool_mat=exp_mat<mat_hpo
      mat_hpo_correctded=mat_hpo+as.data.frame.integer(bool_mat)
      if (enrORdep == "enr") {
        StatMat$'P-value'[hpo] <- stats::fisher.test(mat_hpo_correctded, alternative = "g")$p.value
      } else if (enrORdep == "dep") {
        StatMat$'P-value'[hpo] <- stats::fisher.test(mat_hpo_correctded, alternative = "l")$p.value
      } else if (enrORdep == "both") {
        StatMat$'P-value'[hpo] <- stats::fisher.test(mat_hpo_correctded, alternative = "t")$p.value
      }
      
      # if (enrORdep == "enr") {
      #   StatMat$'P-value'[hpo] <- phyper(n_entered-1,n_world,N_world-n_world,N_entered,lower.tail=F)
      # } else if (enrORdep == "dep") {
      #   StatMat$'P-value'[hpo] <- phyper(n_entered,n_world,N_world-n_world,N_entered,lower.tail=T)
      # } else if (enrORdep == "both") {
      #   if (!is.na(StatMat$Enrichment[hpo])){
      #     if (StatMat$Enrichment[hpo] > 1){
      #       StatMat$'P-value'[hpo] <- phyper(n_entered-1,n_world,N_world-n_world,N_entered,lower.tail=F)
      #     } else {
      #       StatMat$'P-value'[hpo] <- phyper(n_entered,n_world,N_world-n_world,N_entered,lower.tail=T)
      #     }
      #   }
      # 
      # }
    }
  }
  #FDR
  idx <- which(StatMat[,"Observed (n_entered)"] >= minGenes)
  if (length(idx) > 0) {StatMat$FDR[idx] <- p.adjust(StatMat$'P-value'[idx], method = "BH")}
  # idx <- which(StatMat$Enrichment < 1 & StatMat[,"Expected"] >= minGenes)
  # if (length(idx) > 0) {StatMat$FDR[idx] <- p.adjust(StatMat$'P-value'[idx], method = "BH")}
  
  idx <- which(StatMat$FDR < 0.05)
  print(paste("Found ",length(idx)," significant HPO terms",sep=""))
  
  #SAVE OUTPUT
  write.table(StatMat, 
              file=paste(outpath,"HPO_enrich_","v",HPO_dbVer,"_",
                         minGenes,"_",systems,"_","bp_",bodyparts,"_","freq_",
                         TYPorALL,"_",enrORdep,".txt",sep=""), 
              sep="\t", quote = F,na = "", row.names = F, col.names = T)
  
  return(StatMat)
}
