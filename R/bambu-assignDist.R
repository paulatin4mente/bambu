# --- assignReadClasstoTranscripts ---
# Module: Module 4 — Read class to transcript assignment | bambu-assignDist.R
# Called by: bambu.R
# Call count: 1 call, 1 file
#' Create equivilence classes and assign to transcripts
#' @inheritParams bambu
#' @import data.table
#' @noRd
assignReadClasstoTranscripts <- function(readClassList, annotations, isoreParameters, 
                                        verbose, sampleMetadata, demultiplexed,
                                        returnDistTable = FALSE, trackReads = TRUE) {
    if (is.character(readClassList)) readClassList <- readRDS(file = readClassList)
    metadata(readClassList)$readClassDist <- calculateDistTable(readClassList, annotations, isoreParameters, verbose, returnDistTable)
    readClassList <- splitReadClassFiles(readClassList)
    readClassDt <- genEquiRCs(metadata(readClassList)$readClassDist, annotations, verbose) 
    readClassDt$eqClass.match = match(readClassDt$eqClassById,metadata(readClassList)$eqClassById)
    readClassDt <- simplifyNames(readClassDt)
    readClassDt <- readClassDt %>% group_by(eqClassId, gene_sid) %>% 
        mutate(multi_align = length(unique(txid))>1) %>% 
        ungroup() %>% 
        mutate(aval = 1) %>%
        data.table()
    #return non-em counts
    ColData <- generateColData(readClassList, sampleMetadata, demultiplexed)
    quantData <- SummarizedExperiment(assays = SimpleList(
        counts = generateUniqueCounts(readClassDt, metadata(readClassList)$countMatrix, annotations)),
        rowRanges = annotations,
        colData = ColData)
    colnames(quantData) <- ColData$id
    if(sum(metadata(readClassList)$incompatibleCountMatrix)==0){
        metadata(quantData)$incompatibleCounts <- NULL
    }else{
        metadata(quantData)$incompatibleCounts <- generateIncompatibleCounts(metadata(readClassList)$incompatibleCountMatrix, annotations)       
    }
    metadata(quantData)$nonuniqueCounts <- generateNonUniqueCounts(readClassDt, metadata(readClassList)$countMatrix, annotations)
    metadata(quantData)$readClassDt <- readClassDt
    metadata(quantData)$countMatrix <- metadata(readClassList)$countMatrix
    metadata(quantData)$incompatibleCountMatrix <- metadata(readClassList)$incompatibleCountMatrix 
    metadata(quantData)$sampleName <- metadata(readClassList)$sampleData$sampleName 
    if(returnDistTable)
        metadata(quantData)$distTable <- metadata(metadata(readClassList)$readClassDist)$distTableOld

    if(trackReads)
        metadata(quantData)$readToTranscriptMap <- 
            generateReadToTranscriptMap(readClassList, 
                                        metadata(readClassList)$readClassDist, 
                                        annotations)

    return(quantData)     

}

# --- generateUniqueCounts ---
# Module: Module 4 — Read class to transcript assignment | bambu-assignDist.R
# Called by: bambu-assignDist.R
# Call count: 1 call, 1 file
#' Generate unique counts
#' @noRd
generateUniqueCounts <- function(readClassDt, countMatrix, annotations){
    # TODO: [POOR NAMING] x is used to store filtered unique read classes; rename to uniqueReadClassDt or similar
    x <- readClassDt %>% filter(!multi_align & !is.na(eqClass.match))
    uniqueCounts <- countMatrix[x$eqClass.match,]
    uniqueCounts.tx <- sparse.model.matrix(~ factor(x$txid) - 1)
    uniqueCounts <- t(uniqueCounts.tx) %*% uniqueCounts
    rownames(uniqueCounts) <- names(annotations)[match(as.numeric(levels(factor(x$txid))),mcols(annotations)$txid)]
    counts <- sparseMatrix(length(annotations), ncol(uniqueCounts), x = 0)
    rownames(counts) <- names(annotations)
    counts[rownames(uniqueCounts),] <- uniqueCounts
    return(counts)

    # TODO: [UNUSED CODE] the three lines below are unreachable (after return); remove them
    # counts.total = colSums(countMatrix) + colSums(incompatibleCountMatrix)
    # counts.total[counts.total==0] = 1
    # counts.CPM = counts/counts.total * 10^6

}


# --- generateIncompatibleCounts ---
# Module: Module 4 — Read class to transcript assignment | bambu-assignDist.R
# Called by: bambu-assignDist.R
# Call count: 1 call, 1 file
#' Generate incompatible counts
#' @noRd
generateIncompatibleCounts <- function(incompatibleCountMatrix, annotations){
    genes <- levels(factor(unique(mcols(annotations)$GENEID)))
    rownames(incompatibleCountMatrix) <- genes[as.numeric(rownames(incompatibleCountMatrix))]
    geneMat <- sparseMatrix(length(genes), ncol(incompatibleCountMatrix), x = 0)
    rownames(geneMat) <- genes
    geneMat[rownames(incompatibleCountMatrix),] <- incompatibleCountMatrix
    return(geneMat)
}


# --- generateNonUniqueCounts ---
# Module: Module 4 — Read class to transcript assignment | bambu-assignDist.R
# Called by: bambu-assignDist.R
# Call count: 1 call, 1 file
#' Generate non-unique counts
#' @noRd
generateNonUniqueCounts <- function(readClassDt, countMatrix, annotations){
    #fuse multi align RCs by gene
    # TODO: [POOR NAMING] x reused with a different meaning than in generateUniqueCounts (multi-aligned reads); rename to multiAlignReadClassDt
    x <- readClassDt %>% filter(multi_align & !is.na(eqClass.match))
    x <- x %>% distinct(eqClassId, .keep_all = TRUE)
    nonuniqueCounts <- countMatrix[x$eqClass.match,, drop = FALSE]
    if(nrow(x)>1 & length(unique(x$gene_sid))>1){
        nonuniqueCounts.gene <- sparse.model.matrix(~ factor(x$gene_sid) - 1)
        nonuniqueCounts <- t(nonuniqueCounts.gene) %*% nonuniqueCounts
    } else{
        warning("The factor variable 'gene_sid' has only one level. Adjusting output.")
        nonuniqueCounts.gene <- Matrix(1, nrow = nrow(x), ncol = 1, sparse = TRUE)
        nonuniqueCounts <- t(nonuniqueCounts.gene) %*% nonuniqueCounts
    }
    #covert ids into gene ids
    geneids <- as.numeric(levels(factor(x$gene_sid)))
    geneids <- x$txid[match(geneids, x$gene_sid)]
    geneids <- mcols(annotations)$GENEID[as.numeric(geneids)]
    rownames(nonuniqueCounts) <- geneids
    #create matrix for all annotated genes
    genes <- levels(factor(unique(mcols(annotations)$GENEID)))
    geneMat <- sparseMatrix(length(genes), ncol(nonuniqueCounts), x = 0)
    rownames(geneMat) <- genes
    if(!is.null(rownames(nonuniqueCounts))){
      geneMat[rownames(nonuniqueCounts),] <- nonuniqueCounts
    }
    return(geneMat)
}
