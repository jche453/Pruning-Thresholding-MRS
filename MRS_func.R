GenCoMeRegion = function(cmrs = NULL, beta, reference = NULL, Overlap = F){
  if (is.null(cmrs) & is.null(reference)) {
    print("No cmrs file or reference cmrs file. Please input at least one of cmrs file or reference cmrs file")
    exit()
  }
  if (is.null(beta)) {
    print("No DNA methylation (beta) file. Please input DNA methylation (beta) file")
    exit()
  }
  if (Overlap == T){
    if (is.null(cmrs) & !is.null(reference)) {
      print("No cmrs file. Please input cmrs file")
      exit()
    }else if (is.null(reference) & !is.null(cmrs)){
      print("No reference file. Please input reference cmrs file")
      exit()
    }else{
      chmr_refs_flt = unlist(reference,recursive = F) # this is the reference list of 24K CMRs
      cmrs_input_flt = unlist(cmrs,recursive = F) # the evan CMRs
      chmr_refs_prb = unlist(reference) # probes for reference list of 24K CMRs
      cmrs_input_prb = unlist(cmrs) # probes for evan CMRs
      which_ref_ovrlp = which(sapply(chmr_refs_flt, function(x){length(intersect(x,cmrs_input_prb))>0}))
      which_input_NOTovrlp = which(sapply(cmrs_input_flt, function(x){length(intersect(x,chmr_refs_prb))==0}))
      combined_cmr = list()
      combined_cmr = c(combined_cmr, chmr_refs_flt[which_ref_ovrlp], cmrs_input_flt[which_input_NOTovrlp]) 
      cmrsnew = combined_cmr
    }
  }else if (Overlap == F & !is.null(cmrs)){
    print("Use cmrs file to generate CoMeRegion")
    cmrsnew = unlist(cmrs,recursive = F)
  }else if (Overlap == F & !is.null(reference)){
    print("Use reference file to generate CoMeRegion")
    cmrsnew = unlist(reference,recursive = F)
  }
  
  CoMeRegion = NULL 
  for (i in 1: length(cmrsnew)){
    print(paste0("run this ", i, " times"))
    clustercpg = unlist(cmrsnew[i])
    data_cpg = data.frame(clustercpg, CoMethylCluster = i)
    CoMeRegion = rbind(CoMeRegion,data_cpg)
  } 
  print("finished running 1")
  
  CpGnames = colnames(beta)
  print(length(CpGnames))
  
  CpGnames = CpGnames[!(CpGnames %in% CoMeRegion$clustercpg)]
  
  print(length(CpGnames))
  print("finished running 2")
  
  data_cpg = data.frame(clustercpg = CpGnames, 
                        CoMethylCluster = (max(CoMeRegion$CoMethylCluster) +1):(length(CpGnames)  + max(CoMeRegion$CoMethylCluster)))
  CoMeRegion = rbind(CoMeRegion,data_cpg)
  print(dim(CoMeRegion))
  
  return(CoMeRegion)
}

GenMRS = function(beta, SS, Pthred, CoMeRegion, CoMeBack = T, weightSE = F){
  pvalueinfo = matrix(NA, length(Pthred),5)
  colnames(pvalueinfo) = c("Pvalue", "Number of CpG sites", "Numeber of CpG sites after matching with CoMeRegion", " Number of CpG sites after pruning", "Numeber of CpG sites after matching with DNA methylation data")
  if (weightSE == T){
    SS$BETA = SS$BETA/SS$SE
  }
  for (i in Pthred){
    pvalue = 5 * 10 ^ (-i)
    pvalueinfo[i-1,1] =  pvalue
    SS_sub = SS[SS$Pvalue < pvalue, ]
    pvalueinfo[i-1,2] =  nrow(SS_sub)
    
    if (CoMeBack){
      SS_final = merge(SS_sub, CoMeRegion, by.x = "Marker", by.y = "clustercpg")
      SS_final = data.frame(SS_final)
      pvalueinfo[i-1,3] =  nrow(SS_final)
      SS_final = SS_final %>% 
        group_by(CoMethylCluster) %>% 
        slice(which.min(Pvalue))
      pvalueinfo[i-1,4] =  nrow(SS_final)
    }else{
      SS_final = merge(SS_sub, CoMeRegion, by.x = "Marker", by.y = "clustercpg")
      SS_final = data.frame(SS_final)
    }
    betas_final = data.frame(t(beta[,SS_final$Marker, drop = F]))
    betas_final$Marker = rownames(betas_final)
    mat_final = merge(betas_final, SS_final, by = "Marker")
    pvalueinfo[i-1,5] =  nrow(mat_final)
    if (nrow(mat_final) > 0){
      MRS = data.frame(t(as.matrix(mat_final[,2:(ncol(betas_final))])) %*% (mat_final$BETA))
      colnames(MRS) = c(paste0("P",pvalue))
      MRS$ID = rownames(MRS)
      if (i == Pthred[1]){
        MRS_final = MRS
      }else{
        MRS_final = merge(MRS, MRS_final, by = "ID", all = T)
      }
    }
  }
  results = list(pvalueinfo = pvalueinfo, MRS = MRS_final)
  return(results)
}

GenMRS_mediation = function(beta, SS, Pthred){
  pvalueinfo = matrix(NA, length(Pthred),5)
  colnames(pvalueinfo) = c("Pvalue", "Number of CpG sites", "Numeber of CpG sites after matching with CoMeRegion", " Number of CpG sites after pruning", "Numeber of CpG sites after matching with DNA methylation data")
  for (i in Pthred){
    pvalue = 5 * 10 ^ (-i)
    pvalueinfo[i-1,1] =  pvalue
    SS_sub = SS[SS$Pvalue.x < pvalue & SS$Pvalue.y < pvalue, ]
    pvalueinfo[i-1,2] =  nrow(SS_sub)
    
    SS_final = data.frame(SS_sub)
    
    betas_final = data.frame(t(beta[,SS_final$Marker, drop = F]))
    betas_final$Marker = rownames(betas_final)
    mat_final = merge(betas_final, SS_final, by = "Marker")
    pvalueinfo[i-1,5] =  nrow(mat_final)
    if (nrow(mat_final) > 0){
      MRS = data.frame(t(as.matrix(mat_final[,2:(ncol(betas_final))])) %*% (mat_final$BETA))
      colnames(MRS) = c(paste0("P",pvalue))
      MRS$ID = rownames(MRS)
      if (i == Pthred[1]){
        MRS_final = MRS
      }else{
        MRS_final = merge(MRS, MRS_final, by = "ID", all = T)
      }
    }
  }
  results = list(pvalueinfo = pvalueinfo, MRS = MRS_final)
  return(results)
}
