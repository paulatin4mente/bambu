# --- modifyIncompatibleAssignment ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu_utilityFunctions.R
# Call count: 1 call, 1 file
#' modifiy incompatible read classes assignment
#' @import data.table
#' @noRd
modifyIncompatibleAssignment <- function(distTable){
  distTable <- data.table(as.data.frame(distTable))
  distTable[,`:=`(anyCompatible = any(compatible), 
                  anyEqual = any(equal)),
            by = readClassId]
  ## for uncompatible reads, change annotationTxId to gene_id_unidentified
  distTable[which(!anyCompatible), 
            annotationTxId := paste0(GENEID, 
                                     "_unidentified")]
  distTable[,`:=`(anyCompatible = NULL,
                  anyEqual = NULL)]
  distTable <- unique(DataFrame(setDF(distTable)))
  return(distTable)
}


# --- processIncompatibleCounts ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: (not called anywhere)
# Call count: 0 internal calls (exported or not called internally)
#' Process incompatible counts
#' @noRd
processIncompatibleCounts <- function(distTable){
    distTable <- data.table(as.data.frame(distTable))[, 
        .(readClassId, annotationTxId, readCount, GENEID, GENEID.match, GENEID.i, dist,equal)]
    distTable <- distTable[grep("unidentified", annotationTxId)]
    # filter out multiple geneIDs mapped to the same readClass using rowData(se)
    distTable[GENEID.match==TRUE,]
    distTable[, readCount := sum(readCount), by = GENEID]
    counts <- unique(distTable[,.(GENEID, GENEID.i, readCount)])
    setnames(counts, "readCount", "counts")
    return(counts)
}


# --- genEquiRCs ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-assignDist.R
# Call count: 1 call, 1 file
#' This function aims to aggregate RC to equiRCs with the consideration of full
#' or partial alignment status and add empty read class depends on the minimal
#' eqClass from annotations
#' @import data.table
#' @noRd
genEquiRCs <- function(readClassDist, annotations, verbose){
  eqClassCount <- getUniCountPerEquiRC(metadata(readClassDist)$distTable)
  eqClassTable <- addEmptyRC(eqClassCount, annotations)
  # create equiRC id 
  eqClassTable <- eqClassTable %>% 
    group_by(eqClassById) %>%
    mutate(eqClassId = cur_group_id()) %>%
    data.table()

  tx_len <- rbind(data.table(txid = mcols(annotations)$txid,
                             txlen = sum(width(annotations))))
  eqClassTable <- tx_len[eqClassTable, on = "txid"] %>% distinct()
  
  # remove unused columns
  #eqClassTable[, eqClassById := NULL]

  return(eqClassTable)
}

# --- genEquiRCsBasedOnObservedReads ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu_utilityFunctions.R
# Call count: 1 call, 1 file
#' This function formats the distance table obtained from readClass by checking
#' the distance between readClass and transcripts to create equiValence read
#' classes
#' @import data.table
#' @noRd
genEquiRCsBasedOnObservedReads <- function(readClass){
  unlisted_rowranges <- unlist(rowRanges(readClass))
  rcWidth <- data.table(readClassId = rownames(readClass),
                        firstExonWidth =
                          width(unlisted_rowranges[unlisted_rowranges$exon_rank == 1,]),
                        totalWidth = sum(width(rowRanges(readClass))))
  distTable <- data.table(as.data.frame(metadata(readClass)$distTable))[!grepl("unidentified", annotationTxId), .(readClassId, 
                                                                                                                  annotationTxId, readCount, GENEID, dist,equal, compatible, txid)]
  distTable <- rcWidth[distTable, on = "readClassId"]
  # filter out multiple geneIDs mapped to the same readClass using rowData(se)
  compatibleData <- as.data.table(as.data.frame(rowData(readClass)),
                                  keep.rownames = TRUE)
  setnames(compatibleData, old = c("rn", "geneId"),
           new = c("readClassId", "GENEID"))
  distTable <- distTable[compatibleData[ readClassId %in% 
                                           unique(distTable$readClassId), .(readClassId, GENEID)],
                         on = c("readClassId", "GENEID")]
  #here, each transcript should be assigned to one gene only based on isore.estimateDistanceToAnnotation function
  ##this step is very slow, consider to use integers instead of tx_ids
  eqClassByIdList <- createList(distTable$readClassId, distTable$txid*(-1)^distTable$equal)
  distTable[, eqClassById := as.list(eqClassByIdList)]
  return(distTable)
}

# --- createList ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify_utilityFunctions.R
# Call count: 2 calls, 1 file
#' create list
#' @import data.table
#' @noRd
createList <- function(query, subject){
  eqDt <- data.table(query = query, subject = subject)
  eqDt <- eqDt[order(query, subject)]
  eqClassByIdList <- splitAsList(as.integer(eqDt$subject), eqDt$query)
  eqClassByIdList <- unname(eqClassByIdList[as.character(query)])
  return(eqClassByIdList)
}


# --- getUniCountPerEquiRC ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify_utilityFunctions.R
# Call count: 1 call, 1 file
#' get the count and rc_width for each equiRC
#' @import data.table
#' @noRd
getUniCountPerEquiRC <- function(distTable){
  eqClassCount <- distTable %>% 
    group_by(eqClassById) %>%
    mutate(anyEqual = any(equal)) %>%
    select(eqClassById, firstExonWidth,totalWidth, readCount,GENEID,anyEqual) %>% #eqClassByIdTemp,
    distinct() %>%
    mutate(nobs = sum(readCount),
           rcWidth = ifelse(anyEqual, max(totalWidth), 
                            max(firstExonWidth))) %>%
    select(eqClassById,GENEID,nobs,rcWidth) %>% #eqClassByIdTemp,
    ungroup()  %>%
    distinct()
  return(eqClassCount)
}


# --- addEmptyRC ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify_utilityFunctions.R
# Call count: 1 call, 1 file
# add minimal equiRC
#' @import data.table
#' @noRd
addEmptyRC <- function(eqClassCount, annotations){
  minEquiRC <- processMinEquiRC(annotations)
  eqClassCount <- createEqClassToTxMapping(eqClassCount)
  eqClassCountJoin <- full_join(eqClassCount, minEquiRC, by = c("eqClassById","GENEID","txid","equal"))
  eqClassCountJoin[is.na(eqClassCountJoin)] <- 0
  eqClassCount_final <- eqClassCountJoin %>% 
    group_by(eqClassById) %>%
    mutate(nobs = max(nobs),
           rcWidth = max(rcWidth),
           minRC = max(minRC)) %>%
    ungroup() %>%
    distinct()
  return(eqClassCount_final)
}

# --- processMinEquiRC ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify_utilityFunctions.R
# Call count: 1 call, 1 file
# process minEquiRC
#' @import data.table
#' @noRd
processMinEquiRC <- function(annotations){
  minEquiRC <- as.data.frame(mcols(annotations)[,c("eqClassById","GENEID","txid")])
  minEquiRC$eqClassById <- unAsIs(minEquiRC$eqClassById)
  
  minEquiRCTemp <- minEquiRC  %>% 
    mutate(txidTemp = eqClassById) %>% 
    unnest(c(txidTemp)) %>% # split minequirc to txid
    mutate(equal = ifelse(txid == txidTemp, TRUE, FALSE))
  
  eqClassByIdList <- createList(minEquiRCTemp$txid, minEquiRCTemp$txidTemp*(-1)^minEquiRCTemp$equal)
  minEquiRCTemp$eqClassById <- as.list(eqClassByIdList)
  
  minEquiRCTemp <- minEquiRCTemp %>%
    mutate(txid = txidTemp, eqClassByIdTemp = NULL, txidTemp = NULL) %>%
    distinct()
  
  minEquiRC <- minEquiRC %>% 
    mutate(txidTemp = eqClassById) %>% 
    unnest(c(txidTemp)) %>%
    mutate(minRC = 1, equal = FALSE,  txid = txidTemp, txidTemp = NULL)
  
  minEquiRC <- bind_rows(minEquiRC, minEquiRCTemp) 
  return(minEquiRC)
}

# --- unAsIs ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify_utilityFunctions.R
# Call count: 1 call, 1 file
#' Function to get rid of AsIs class so that group_by can be used on eqClassById column in mcols(annotations)
#' credit to https://stackoverflow.com/questions/12865218/getting-rid-of-asis-class-attribute
#' @noRd
unAsIs <- function(X) {
  if("AsIs" %in% class(X)) {
    class(X) <- class(X)[-match("AsIs", class(X))]
  }
  X
}


# --- createEqClassToTxMapping ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify_utilityFunctions.R
# Call count: 1 call, 1 file
#' Create eqClass to tx mapping based on eqClassById
#' @import tidyr
#' @noRd
createEqClassToTxMapping <- function(eqClassTable){
  eqClassTable_unnest <- eqClassTable %>% 
    mutate(txid = eqClassById) %>% 
    unnest(c(txid)) %>%
    mutate(equal = ifelse(txid < 0,TRUE,FALSE)) %>%
    mutate(txid = abs(txid))
  return(eqClassTable_unnest)
}

# --- addAval ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify.R
# Call count: 1 call, 1 file
#' Add A matrix for total, full-length, unique
#' @noRd
addAval <- function(readClassDt, emParameters, verbose){
#   if (is.null(readClassDt)) {
#     stop("Input object is missing.")
#   } else if (any(!(c("GENEID", "txid", "eqClassId","nobs") %in% 
#                    colnames(readClassDt)))) {
#     stop("Columns GENEID, txid, eqClassId, nobs,
#             are missing from object.")
#   }
  ## ----step2: match to simple numbers to increase claculation efficiency
  #readClassDt <- simplifyNames(readClassDt)
  d_mode <- emParameters[["degradationBias"]]
  start.ptm <- proc.time()
  if (d_mode) {
    d_rateOut <- calculateDegradationRate(readClassDt)
    readClassDt <- modifyAvaluewithDegradation_rate(readClassDt, 
                                                    d_rateOut[1], d_mode = d_mode)
  }else{
    d_rateOut <- rep(NA,2)
    readClassDt$aval = 1
  }
  end.ptm <- proc.time()
  if (verbose) message("Finished estimate degradation bias in ",
                       round((end.ptm - start.ptm)[3] / 60, 3), " mins.")
  removeList <- removeUnObservedGenes(readClassDt)
  readClassDt <- removeList[[1]] # keep only observed genes for estimation
  outList <- removeList[[2]] #for unobserved genes, set estimates to 0 
  readClassDt_withGeneCount <- select(readClassDt, gene_sid, eqClassId, nobs) %>%
    unique() %>%
    group_by(gene_sid) %>%
    mutate(K = sum(nobs), n.obs=nobs/K, nobs = NULL) %>% ## check if this is unique by eqClassId
    ungroup() %>%
    distinct() %>%
    right_join(readClassDt, by = c("gene_sid","eqClassId"))
#   readClassDt_withGeneCount <- readClassDt %>%
#     group_by(gene_sid, eqClassId, nobs) %>% mutate(n = n(), n = replace(n,1,1)) %>%
#     ungroup() %>% group_by(gene_sid) %>%
#     mutate(K = sum(nobs), n.obs=(nobs)/(K)) %>%
#     ungroup()
  return(list(readClassDt_withGeneCount,outList))
}

# --- simplifyNames ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-assignDist.R
# Call count: 1 call, 1 file
#' This function converts transcript, gene, and read class names to simple
#' integers for more efficient computation
#' @import data.table
#' @noRd
simplifyNames <- function(readClassDt){
  readClassDt <- as.data.table(readClassDt)
  readClassDt[, gene_sid := match(GENEID, unique(readClassDt$GENEID))]
  readClassDt[, `:=`(GENEID = NULL)]
  return(readClassDt)
}


# --- calculateDegradationRate ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify_utilityFunctions.R
# Call count: 1 call, 1 file
#' Calculate degradation rate based on equiRC read counts
#' @import data.table
#' @noRd
calculateDegradationRate <- function(readClassDt){
  rcCount <- unique(readClassDt[, .(gene_sid,eqClassId, nobs)])
  rcCountPar <-
    unique(readClassDt[which(!equal), .(gene_sid,eqClassId, nobs)])
  geneCount <- unique(rcCount[, list(nobs = sum(nobs)), by = gene_sid])
  geneCountPar <- unique(rcCountPar[, list(dObs = sum(nobs)), by = gene_sid])
  txLength <- unique(readClassDt[, .(gene_sid, txid, txlen)])
  geneLength <- 
    unique(txLength[, list(gene_len = max(txlen)), by = gene_sid])
  geneCountLength <- unique(geneLength[geneCount, on = "gene_sid"])
  geneCountLength <- unique(geneCountPar[geneCountLength, on = "gene_sid"])
  geneCountLength[, d_rate := dObs/nobs]
  if (length(which(geneCountLength$nobs >= 30 & 
                   ((geneCountLength$nobs - geneCountLength$dObs) >= 5))) == 0) {
    # message("There is not enough read count and full length coverage!
    #         Hence degradation rate is estimated using all data!")
  } else {
    geneCountLength <- geneCountLength[nobs >= 30 & ((nobs - dObs) >= 5)]
  }
  d_rate <- median(geneCountLength$d_rate * 1000/geneCountLength$gene_len,
                   na.rm = TRUE)
  return(c(d_rate, nrow(geneCountLength)))
}



# --- modifyAvaluewithDegradation_rate ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify_utilityFunctions.R
# Call count: 1 call, 1 file
#' This function generates a_mat values for all transcripts
#' @import data.table
#' @noRd
modifyAvaluewithDegradation_rate <- function(tmp, d_rate, d_mode){
  tmp[, multi_align := (length(unique(txid)) > 1),
      by = list(eqClassId, gene_sid)]
  if (!d_mode) {
    tmp[, aval := 1]
    return(tmp)
  }
  tmp[which(multi_align) , aval := ifelse(equal, 1 -
                                            sum(.SD[which(!equal)]$rcWidth*d_rate/1000),
                                          rcWidth*d_rate/1000), by = list(gene_sid,txid)]
  if(is.na(d_rate)) d_rate = 0
  if (d_rate == 0) {
    tmp[, par_status := all(!equal & multi_align),
        by = list(eqClassId, gene_sid)]
    tmp[which(par_status), aval := 0.01]
  }
  tmp[, aval := pmax(pmin(aval,1),0)] #d_rate should be contained to 0-1
  tmp[multi_align & equal,
      aval := pmin(1,pmax(aval,rcWidth*d_rate/1000))]
  tmp[which(!multi_align), aval := 1]
  return(tmp)
}


# --- removeUnObservedGenes ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify_utilityFunctions.R
# Call count: 1 call, 1 file
#' @import data.table
#' @noRd
removeUnObservedGenes <- function(readClassDt){
    uoGenes <- unique(readClassDt[,.I[sum(nobs) == 0], by = gene_sid]$gene_sid)
  if (length(uoGenes) > 0) {
    uo_txGeneDt <- 
      unique(readClassDt[(gene_sid %in% uoGenes),.(txid,gene_sid)])
    readClassDt <- readClassDt[!(gene_sid %in% uoGenes)]
    outList <- data.table(txid = uo_txGeneDt$txid,
                          counts = 0, 
                          fullLengthCounts = 0,
                          uniqueCounts = 0)
  }else{
    outList <- NULL
  }
  return(list(readClassDt, outList))
}



# --- initialiseOutput ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify.R
# Call count: 1 call, 1 file
#' This function initialises the final estimates with default values
#' @import data.table
#' @noRd
initialiseOutput <- function(readClassDt){
  return(unique(data.table(txid = readClassDt$txid,
                           counts = 0,
                           fullLengthCounts = 0,
                           uniqueCounts = 0),by = NULL))
}


# --- filterTxRc ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify.R
# Call count: 1 call, 1 file
#' filter transcripts without read support
#' @noRd
filterTxRc <- function(readClassDt){
  readClassDt <- group_by(readClassDt,gene_sid, txid) %>% 
    filter(sum(nobs)>0) %>% ungroup() %>%
    arrange(gene_sid, txid, eqClassId) %>%
    mutate(uniqueAval =  aval * (!multi_align),
           fullAval = aval * equal,
           multi_align = NULL,
           equal = NULL) %>%
    data.table()
  return(readClassDt)
}


# --- assignGroups ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify.R
# Call count: 1 call, 1 file
#' Assign internal groups for grouped fast processing
#' @noRd
assignGroups <- function(readClassDt){
  # further devide genes into groups to improve process efficiency
  readClassDt[, tr_dimension := length(unique(txid))*length(unique(eqClassId)), by = gene_sid]
  trDimensionDt <- unique(readClassDt[,.(gene_sid, tr_dimension)])
  trDimensionDt[order(tr_dimension), cumN := cumsum(tr_dimension)]
  trDimensionDt[, gene_grp_id := cumN %/% 1000 + 1]
  readClassDt <- trDimensionDt[readClassDt, on = c("gene_sid","tr_dimension")]
  readClassDt[, `:=`(cumN = NULL, tr_dimension = NULL)]
  return(readClassDt)
}
# --- getInputList ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify.R
# Call count: 1 call, 1 file
#'
#' @noRd
getInputList <- function(readClassDt){
  nObsVec <- unique(readClassDt[,.(gene_grp_id, gene_sid, eqClassId, K, n.obs)])[order(gene_grp_id, gene_sid, eqClassId)]
  nObsList <- splitAsList(nObsVec$n.obs,nObsVec$gene_grp_id)
  inputRcDt <- unique(readClassDt[,.(gene_grp_id)][order(gene_grp_id)])
  inputRcDt[, nObs_list := as.list(unname(nObsList))]
  KList <- splitAsList(nObsVec$K,nObsVec$gene_grp_id)
  inputRcDt[, K_list := as.list(unname(KList))]
  txidsVec <- unique(readClassDt[,.(gene_grp_id, txid)])[order(gene_grp_id,txid)]
  txidsList <- splitAsList(txidsVec$txid,txidsVec$gene_grp_id)
  inputRcDt[, txids_list := as.list(unname(txidsList))]
  inputRcDt <- split(inputRcDt, by = "gene_grp_id")
  return(inputRcDt)
}


# --- abundance_quantification ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify.R
# Call count: 1 call, 1 file
#' Nanopore transcript abundance quantification
#' @title transcript_abundance_quantification
#' @param readClassDt A \code{data.table} with columns
#' @importFrom BiocParallel bpparam bplapply
#' @noRd
abundance_quantification <- function(inputRcDt, readClassDt,
                                     maxiter = 20000, conv = 10^(-2), minvalue = 10^(-8)) {
        emResultsList <- lapply(as.list(names(inputRcDt)),
                                run_parallel,
                                conv = conv,
                                minvalue = minvalue, 
                                maxiter = maxiter,
                                inputRcDt = inputRcDt,
                                readClassDt = readClassDt
        )
    estimates <- do.call("rbind", emResultsList)
    return(estimates)
}


# --- run_parallel ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: (not called anywhere)
# Call count: 0 internal calls (exported or not called internally)
#' function to run in parallel
#' @title run_parallel
#' @param g the serial id of gene
#' @inheritParams abundance_quantification
#' @importFrom methods is
#' @import data.table
#' @noRd
run_parallel <- function(g, conv, minvalue, maxiter, inputRcDt, readClassDt) {
    input_g <- inputRcDt[[g]]  
    K <- input_g$K_list[[1]]
    n.obs <- input_g$nObs_list[[1]]
    txids <- input_g$txids_list[[1]]
    rcMat <- readClassDt[[g]]
    total_mat <- getAMat(rcMat, by = "aval")
    full_mat <- getAMat(rcMat, by = "fullAval")
    unique_mat <- getAMat(rcMat, by = "uniqueAval")
     if (is(total_mat,"numeric")&(nrow(total_mat)==1)) {
        outEst <- cbind(sum(K * n.obs * total_mat),
                       sum(K * n.obs * full_mat),
                       sum(K * n.obs * unique_mat),
                       unname(txids))
    }else{
        est_output <- emWithL1(A = total_mat, A_full = full_mat, A_unique = unique_mat, Y = n.obs, K = K,
                                maxiter = maxiter,
                               minvalue = minvalue, conv = conv)
        outEst <- cbind(t(est_output[["theta"]]),unname(txids))
    }
    return(outEst)
}


# --- getAMat ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify_utilityFunctions.R
# Call count: 3 calls, 1 file
#' Get the A matrix
#' @noRd
getAMat <- function(rcMat, by = "aval"){
    a_mat <- dcast(rcMat, gene_sid + eqClassId ~ txid, value.var = by, fill = 0)
    a_mat[, `:=`(gene_sid = NULL, eqClassId = NULL)]
    return(t(as.matrix(a_mat)))
}



# --- modifyQuantOut ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify.R
# Call count: 1 call, 1 file
#' Modify default quant output using estimated outputs
#' @import data.table
#' @noRd
modifyQuantOut <- function(outEst, outIni){
    outIni[match(as.integer(outEst[,4]), txid),counts := outEst[,1]]
    outIni[match(as.integer(outEst[,4]), txid), fullLengthCounts := outEst[,2]]
    outIni[match(as.integer(outEst[,4]), txid), uniqueCounts := outEst[,3]]
    return(outIni)
}


# --- removeDuplicates ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify.R
# Call count: 1 call, 1 file
#' Remove duplicate transcript counts originated from multiple genes
#' @import data.table
#' @noRd
removeDuplicates <- function(counts){
    counts_final <- unique(counts[, list(counts = sum(counts),
                                         fullLengthCounts = sum(fullLengthCounts),
                                         uniqueCounts = sum(uniqueCounts)), 
                                  by = txid],by = NULL)
    return(counts_final)
}





# --- generateReadToTranscriptMap ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-assignDist.R
# Call count: 1 call, 1 file
#' Generate read to transcript mapping
#' @noRd
generateReadToTranscriptMap <- function(readClass, distTable, annotations){
    if(!is.null(metadata(readClass)$readNames)) { 
        read_id <- metadata(readClass)$readNames
    } else { 
        read_id <- metadata(readClass)$readId
    }
  #unpack and reverse the read class to read id relationship
  readOrder <- order(unlist(rowData(readClass)$readIds))
  lens <- lengths(rowData(readClass)$readIds)
  rcIndex <- seq_along(readClass)
  readToRC <- rep(rcIndex, lens)[readOrder]
  read_id <- read_id[(match(unlist(rowData(readClass)$readIds)[readOrder], 
                           metadata(readClass)$readId))]
  #get annotation indexs
  distTable <- metadata(distTable)$distTable
  distTable$annotationTxId <- match(distTable$annotationTxId, names(annotations))
  #match read classes with transcripts
  readClass_id <- rownames(readClass)[readToRC]
  distTable$exClassById <- NULL
  equalMatches <- as_tibble(distTable) %>%
    filter(equal) %>% 
    group_by(readClassId) %>% 
      summarise(annotationTxIds = list(annotationTxId))
  equalMatches <- equalMatches$annotationTxIds[match(readClass_id,
                        equalMatches$readClassId)]
  compatibleMatches <- as_tibble(distTable) %>% 
    filter(!equal & compatible) %>% 
    group_by(readClassId) %>% summarise(annotationTxIds = list(annotationTxId))
  compatibleMatches <- compatibleMatches$annotationTxIds[match(readClass_id,
                            compatibleMatches$readClassId)]
  readToTranscriptMap <- tibble(readId=read_id, equalMatches = equalMatches, 
                             compatibleMatches = compatibleMatches)
  return(readToTranscriptMap)
}

# --- calculateEqClassCounts ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: (not called anywhere)
# Call count: 0 internal calls (exported or not called internally)
#' Get counts of equivilent classes from a distTable and match to a readClassDt
#' @noRd
calculateEqClassCounts <- function(distTable, readClassDt){
        eqClasses <- distTable %>% group_by(eqClassById) %>%
                         mutate(anyEqual = any(equal)) %>%
                         select(eqClassById, firstExonWidth,totalWidth,
                                readCount,GENEID,anyEqual) %>% #eqClassByIdTemp,
                         distinct() %>%
                         mutate(nobs = sum(readCount),
                             rcWidth = ifelse(anyEqual, max(totalWidth), 
                             max(firstExonWidth))) %>%
                         select(eqClassById,GENEID,nobs,rcWidth) %>% 
                         ungroup() %>% distinct()
        eqCounts <- eqClasses$nobs[match(readClassDt$eqClassById,eqClasses$eqClassById)]
        eqCounts[is.na(eqCounts)] <- 0
        return(eqCounts)
}

# --- calculateCPM ---
# Module: Module 4 — Read class to transcript assignment | bambu-quantify_utilityFunctions.R
# Called by: bambu-quantify.R
# Call count: 1 call, 1 file
#' calculate CPM post estimation
#' @noRd
calculateCPM <- function(compatibleCounts, incompatibleCounts){
    totalCount <- sum(compatibleCounts$counts)+sum(incompatibleCounts$counts)
    compatibleCounts[, `:=`(CPM = counts / totalCount * (10^6))]
    return(compatibleCounts)
}

#' @useDynLib bambu, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


#' @noRd
.onUnload <- function(libpath) {
    library.dynam.unload("bambu", libpath)
}

