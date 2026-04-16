# --- bambu.processReads ---
# Module: Module 2 — Read processing (per sample) | bambu-processReads.R
# Called by: bambu.R
# Call count: 1 calls, 1 files
#' process reads
#' @param reads path to BAM file(s)
#' @param annotations path to GTF file or TxDb object
#' @param genomeSequence path to FA file or BSgenome object
#' @param readClass.outputDir path to readClass output directory
#' @param yieldSize yieldSize
#' @param bpParameters BioParallel parameter
#' @param stranded stranded
#' @param verbose verbose
#' @importFrom Rsamtools yieldSize BamFileList yieldSize<-
#' @importFrom methods is
#' @importFrom BiocParallel bplapply
#' @importFrom BiocGenerics basename
#' @noRd
bambu.processReads <- function(reads, annotations, genomeSequence,
    readClass.outputDir=NULL, yieldSize=1000000, bpParameters, 
    stranded=FALSE, verbose=FALSE, isoreParameters = setIsoreParameters(NULL),
    processByChromosome = FALSE, processByBam = TRUE, trackReads = trackReads, fusionMode = fusionMode, 
    demultiplexed = FALSE, cleanReads = FALSE, dedupUMI = FALSE, sampleNames = NULL, barcodesToFilter = NULL) {
    genomeSequence <- checkInputSequence(genomeSequence) #TODO (JG) [validate-input]  move to bambu() input validation
    # ===# create BamFileList object from character #===#
    if (is(reads, "BamFile")) { #TODO (JG) [validate-input] should be done in bambu in validate input section, minimisa arguments for this call, fix type to BamFileList with names
        if (!is.null(yieldSize)) {
            yieldSize(reads) <- yieldSize
        } else {
            yieldSize <- yieldSize(reads)
        }
        reads <- BamFileList(reads)
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(reads))
    } else if (is(reads, "BamFileList")) {
        if (!is.null(yieldSize)) {
            yieldSize(reads) <- yieldSize
        } else {
            yieldSize <- min(yieldSize(reads))
        }
    } else if (any(!grepl("\\.bam$", reads))) {
        stop("Bam file is missing from arguments.")#TODO (JG) [validate-input] should be done in bambu in validate input section
    } else {
        if (is.null(yieldSize)) yieldSize <- NA
        reads <- BamFileList(reads, yieldSize = yieldSize)
        names(reads) <- tools::file_path_sans_ext(BiocGenerics::basename(reads))
    }
    if(!is.null(sampleNames)){ #TODO (JG) [validate-input] should be done in bambu in validate input section
        # TODO: [BUG] operator precedence error: length(sampleNames==length(reads)) compares
        # sampleNames to an integer first, producing a logical vector whose length() is always
        # >0, so this condition is always TRUE. Should be: length(sampleNames)==length(reads)
        if(length(sampleNames==length(reads))){#TODO (JG) [bug] see above
            names(reads) <- sampleNames
        } else{
            message("Not enough provided sample names. Using them in order of inputted files and the remaining files will use the file names")
            names(reads)[seq_along(sampleNames)] <- sampleNames
        }
    }
    min.readCount <- isoreParameters[["min.readCount"]]
    fitReadClassModel <- isoreParameters[["fitReadClassModel"]]
    defaultModels <- isoreParameters[["defaultModels"]]
    returnModel <- isoreParameters[["returnModel"]]
    min.exonOverlap <- isoreParameters[["min.exonOverlap"]]

    if(processByBam){ #TODO (JG) [rewrite-processByBam] processByBam can be default to TRUE, possibly remove the part below to combine read classes across files. redundant, difficult to maintain?
        readClassList <- bplapply(seq_along(reads), function(i) { #TODO (JG) [rewrite-processByBam] index is hardcoded to 1 and not used here, and also not used in function below
            bambu.processReadsByFile(bam.file = reads[i],
            genomeSequence = genomeSequence,annotations = annotations,
            stranded = stranded, min.readCount = min.readCount, 
            fitReadClassModel = fitReadClassModel, min.exonOverlap = min.exonOverlap, 
            defaultModels = defaultModels, returnModel = returnModel, verbose = verbose, 
            processByChromosome = processByChromosome, trackReads = trackReads, fusionMode = fusionMode, 
            demultiplexed = demultiplexed, cleanReads = cleanReads, dedupUMI = dedupUMI, index = 1, barcodesToFilter = barcodesToFilter)},
            BPPARAM = bpParameters)
    } else {
        readGrgList <- bplapply(seq_along(reads), function(i) {
            bambu.readsByFile(bam.file = reads[i],
            genomeSequence = genomeSequence,annotations = annotations,
            stranded = stranded, min.readCount = min.readCount, 
            fitReadClassModel = fitReadClassModel, min.exonOverlap = min.exonOverlap, 
            defaultModels = defaultModels, returnModel = returnModel, verbose = verbose, 
            trackReads = trackReads, fusionMode = fusionMode, 
            demultiplexed = demultiplexed, cleanReads = cleanReads, dedupUMI = dedupUMI, index = i, barcodesToFilter = barcodesToFilter)},
            BPPARAM = bpParameters)
        sampleNames <- as.numeric(as.factor(sampleNames))
        for(i in seq_along(readGrgList)){
            if(!isFALSE(demultiplexed)){
                mcols(readGrgList[[i]])$CB <- paste0(names(reads)[i], '_', mcols(readGrgList[[i]])$CB)
            } else{
                mcols(readGrgList[[i]])$CB <- sampleNames[i]
            }
            
            mcols(readGrgList[[i]])$CB <- as.factor(mcols(readGrgList[[i]])$CB)
            
        }
        readGrgList <- do.call(c, readGrgList)    
        mcols(readGrgList)$id <- seq_along(readGrgList) 
        if(!isFALSE(demultiplexed)){ 
          mcols(readGrgList)$sampleID <- as.numeric(mcols(readGrgList)$CB)
        } else {
          mcols(readGrgList)$sampleID <- i
        }
        readClassList <- constructReadClasses(readGrgList, genomeSequence = genomeSequence,annotations = annotations,
            stranded = stranded, min.readCount = min.readCount, 
            fitReadClassModel = fitReadClassModel, min.exonOverlap = min.exonOverlap, 
            defaultModels = defaultModels, returnModel = returnModel, verbose = verbose, 
            processByChromosome = processByChromosome, trackReads = trackReads, fusionMode = fusionMode)
        metadata(readClassList)$samples <- names(reads)
        metadata(readClassList)$sampleNames <- names(reads)
        if(!isFALSE(demultiplexed)) metadata(readClassList)$samples <- levels(mcols(readGrgList)$CB)
        readClassList <- list(readClassList)
    }
        
    if (!is.null(readClass.outputDir)) {
        for(i in seq_along(readClassList)){
            readClassFile <- "combinedSamples"
            readClassFile <- BiocFileCache::bfcnew(BiocFileCache::BiocFileCache(
                readClass.outputDir, ask = FALSE),
                paste0(readClassFile,"_readClassSe"), ext = ".rds")
            saveRDS(readClassList[[i]], file = readClassFile)
            readClassList[[i]] <- readClassFile
        }
    }
    #TODO don't output list, current there because discovery needs it
    return(readClassList)
}

# --- bambu.processReadsByFile ---
# Module: Module 2 — Read processing (per sample) | bambu-processReads.R
# Called by: bambu-processReads.R
# Call count: 1 calls, 1 files
#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @importFrom GenomeInfoDb seqlevels seqlevels<- keepSeqlevels
#' @noRd
bambu.processReadsByFile <- function(bam.file, genomeSequence, annotations,
    yieldSize = NULL, stranded = FALSE, min.readCount = 2, 
    fitReadClassModel = TRUE, min.exonOverlap = 10, defaultModels = NULL, returnModel = FALSE, 
    verbose = FALSE, processByChromosome = FALSE, trackReads = FALSE, fusionMode = FALSE, demultiplexed = FALSE, 
    cleanReads = FALSE, dedupUMI = FALSE, index = 0, barcodesToFilter = NULL) {
    if(verbose) message(names(bam.file)[1])
    readGrgList <- prepareDataFromBam(bam.file[[1]], verbose = verbose, yieldSize = yieldSize, use.names = trackReads, demultiplexed = demultiplexed, cleanReads = cleanReads, dedupUMI = dedupUMI)
    if(verbose) message(paste0("Number of alignments/reads: ",length(readGrgList)))
    warnings <- c() #TODO (JG) [warnings] need to be implemented
    if(!is.null(barcodesToFilter) & !isFALSE(demultiplexed))
        readGrgList <- readGrgList[!(mcols(readGrgList)$CB %in% barcodesToFilter)]
    warnings <- seqlevelCheckReadsAnnotation(readGrgList, annotations)
    if(verbose & length(warnings) > 0) warning(paste(warnings,collapse = "\n"))
    #check seqlevels for consistency, drop ranges not present in genomeSequence
    refSeqLevels <- seqlevels(genomeSequence) 
    if (!all(seqlevels(readGrgList) %in% refSeqLevels)) {
        refSeqLevels <- intersect(refSeqLevels, seqlevels(readGrgList))
        if (!all(seqlevels(annotations) %in% refSeqLevels)&(!(length(annotations)==0))) {
            refSeqLevels <- intersect(refSeqLevels, seqlevels(annotations))
            warningText <- paste0("not all chromosomes from annotations present in ", 
            "reference genome sequence, annotations without reference genomic sequence ",
            "are dropped")
            warnings <- c(warnings, warningText)
            if(verbose) warning(warningText)
            annotations <- keepSeqlevels(annotations, value = refSeqLevels,
                                         pruning.mode = "coarse")
        }
        warningText <- paste0("not all chromosomes from reads present in reference ",
        "genome sequence, reads without reference chromosome sequence are dropped")
        warnings <- c(warnings, warningText)
        if(verbose) warning(warningText)
        readGrgList <- keepSeqlevels(readGrgList, value =  refSeqLevels,
                                     pruning.mode = "coarse")
        # reassign Ids after seqlevels are dropped
        mcols(readGrgList)$id <- seq_along(readGrgList) #TODO (JG) [unused-code] this line is redundant with the line below
    }
    #removes reads that are outside genome coordinates
    badReads <- which(max(end(ranges(readGrgList)))>
                         seqlengths(genomeSequence)[as.character(getChrFromGrList(readGrgList))])
    if(length(badReads) > 0 ){
        readGrgList <- readGrgList[-badReads]
        warningText <- paste0(length(badReads), " reads are mapped outside the provided ",
                       "genomic regions. These reads will be dropped. Check you are using the ",
                       "same genome used for the alignment")
        warnings <- c(warnings, warningText)
        if(verbose) warning(warningText)
    }
    if(length(readGrgList) == 0)
        stop("No reads left after filtering.")

    mcols(readGrgList)$id <- seq_along(readGrgList) 

    if(!isFALSE(demultiplexed)){ 
        mcols(readGrgList)$sampleID <- as.numeric(mcols(readGrgList)$CB)
    } else {
        mcols(readGrgList)$sampleID <- index #TODO (JG) [rewrite-processByBam]index option can be removed if t seems to be hardcoded to 1, as this option can't be changed?
    }
        
    # construct read classes for each chromosome seperately
    if(processByChromosome){
        se <- lowMemoryConstructReadClasses(readGrgList, genomeSequence,
                                                      annotations, stranded, verbose,bam.file)
    } else{
        unlisted_junctions <- unlistIntrons(readGrgList, use.ids = TRUE)
        uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions,
                                                         annotations,genomeSequence, stranded = stranded, verbose = verbose)
        # TODO: [OTHER] runName = "TODO" is a placeholder; replace with a meaningful sample/run identifier
        se <- isore.constructReadClasses(readGrgList,
                                              unlisted_junctions, uniqueJunctions, runName = "TODO",
                                              annotations, stranded, verbose)

    }

    metadata(se)$warnings <- warnings
    if(trackReads){
        metadata(se)$readNames <- names(readGrgList)
        metadata(se)$readId <- mcols(readGrgList)$id
    }
    refSeqLevels <- seqlevels(genomeSequence)
    GenomeInfoDb::seqlevels(se) <- refSeqLevels
    # create SE object with reconstructed readClasses
    se <- scoreReadClasses(se, genomeSequence, annotations, 
                             defaultModels = defaultModels,
                             fit = fitReadClassModel,
                             returnModel = returnModel,
                             min.readCount = min.readCount,
                             min.exonOverlap = min.exonOverlap,
                             fusionMode = fusionMode,
                             verbose = verbose)

    if (demultiplexed) {
        barcodes <- levels(mcols(readGrgList)$CB)
        metadata(se)$sampleData <- tibble(
          id = paste(names(bam.file)[1], barcodes, sep = '_'),
          sampleName = names(bam.file)[1],
          barcode = barcodes
        )
    } else{
        metadata(se)$sampleData <- tibble(
          id = names(bam.file)[1],
          sampleName = names(bam.file)[1]
        )
    }

    return(se)
}

# --- bambu.readsByFile ---
# Module: Module 2 — Read processing (per sample) | bambu-processReads.R
# Called by: bambu-processReads.R
# Call count: 1 calls, 1 files
#' Preprocess bam files and save read class files
#' @inheritParams bambu
#' @importFrom GenomeInfoDb seqlevels seqlevels<- keepSeqlevels
#' @noRd
bambu.readsByFile <- function(bam.file, genomeSequence, annotations,
    yieldSize = NULL, stranded = FALSE, min.readCount = 2, 
    fitReadClassModel = TRUE, min.exonOverlap = 10, defaultModels = NULL, returnModel = FALSE, 
    verbose = FALSE, trackReads = FALSE, fusionMode = FALSE, demultiplexed = FALSE, 
    cleanReads = TRUE, dedupUMI = FALSE, index = 0, barcodesToFilter = NULL) {
    readGrgList <- prepareDataFromBam(bam.file[[1]], verbose = verbose, yieldSize = yieldSize, use.names = trackReads, demultiplexed = demultiplexed, cleanReads = cleanReads, dedupUMI = dedupUMI)
    
    if(!is.null(barcodesToFilter) & !isFALSE(demultiplexed)) readGrgList <- readGrgList[!mcols(readGrgList)$CB %in% barcodesToFilter]
    
    if(verbose) message("Number of alignments/reads: ",length(readGrgList))
    
    warnings <- c()
    warnings <- seqlevelCheckReadsAnnotation(readGrgList, annotations)
    
    if(verbose & length(warnings) > 0) warning(paste(warnings,collapse = "\n"))
    #check seqlevels for consistency, drop ranges not present in genomeSequence
    refSeqLevels <- seqlevels(genomeSequence)
    if (!all(seqlevels(readGrgList) %in% refSeqLevels)) {
        refSeqLevels <- intersect(refSeqLevels, seqlevels(readGrgList))
        if (!all(seqlevels(annotations) %in% refSeqLevels)&(!(length(annotations)==0))) {
          refSeqLevels <- intersect(refSeqLevels, seqlevels(annotations))
          warningText <- paste0("not all chromosomes from annotations present in ", 
                               "reference genome sequence, annotations without reference genomic sequence ",
                               "are dropped")
          warnings <- c(warnings, warningText)
          if(verbose) warning(warningText)
          annotations <- keepSeqlevels(annotations, value = refSeqLevels,
                                       pruning.mode = "coarse")
        }
        warningText <- paste0("not all chromosomes from reads present in reference ",
                             "genome sequence, reads without reference chromosome sequence are dropped")
        warnings <- c(warnings, warningText)
        if(verbose) warning(warningText)
        readGrgList <- keepSeqlevels(readGrgList, value =  refSeqLevels,
                                     pruning.mode = "coarse")
        # reassign Ids after seqlevels are dropped
        mcols(readGrgList)$id <- seq_along(readGrgList) 
    }
    #removes reads that are outside genome coordinates
    badReads <- which(max(end(ranges(readGrgList)))>=
                         seqlengths(genomeSequence)[as.character(getChrFromGrList(readGrgList))])
      if(length(badReads) > 0 ){
        readGrgList <- readGrgList[-badReads]
        warningText <- paste0(length(badReads), " reads are mapped outside the provided ",
                             "genomic regions. These reads will be dropped. Check you are using the ",
                             "same genome used for the alignment")
        warnings <- c(warnings, warningText)
        if(verbose) warning(warningText)
      }
      
      ### add ### 
      # reassign Ids after seqlevels are dropped
      mcols(readGrgList)$id <- seq_along(readGrgList) 
      ### add ###
      if(verbose) message("Number of post-filter alignments/reads: ",length(readGrgList))
      if(length(readGrgList) == 0)
        stop("No reads left after filtering.")
      
      ## add ###
      #if (isTRUE(demultiplexed)){
      #  cellBarcodeAssign <- tibble(index = mcols(readGrgList)$id, CB = mcols(readGrgList)$CB) %>% nest(.by = "CB")

        # if (!dir.exists("CB")){
        #   dir.create("CB")
        # } else{
        #   unlink(paste("CB", "*", sep = "/"))
        # }
        
        # invisible(lapply(seq(nrow(cellBarcodeAssign)),
        #           function(x){saveRDS(readGrgList[pull(cellBarcodeAssign$data[[x]])], paste0("CB/", cellBarcodeAssign$CB[[x]],".rds"))}))
      #} 
    return(readGrgList)
}

# --- constructReadClasses ---
# Module: Module 2 — Read processing (per sample) | bambu-processReads.R
# Called by: bambu-processReads.R
# Call count: 4 calls, 1 files
#' Construct read classes
#' @noRd
constructReadClasses <- function(readGrgList, genomeSequence, annotations,
    stranded = FALSE, min.readCount = 2, 
    fitReadClassModel = TRUE, min.exonOverlap = 10, defaultModels = NULL, returnModel = FALSE, 
    verbose = FALSE, processByChromosome = FALSE, trackReads = FALSE, fusionMode = FALSE){
    
    if(processByChromosome){
        # construct read classes for each chromosome seperately
        # TODO: [OTHER] "TODO" passed as runName is a placeholder; replace with a meaningful sample/run identifier
        se <- lowMemoryConstructReadClasses(readGrgList, genomeSequence,
                                            annotations, stranded, verbose,"TODO", fusionMode)
    } else{
        unlisted_junctions <- unlistIntrons(readGrgList, use.ids = TRUE)
        uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions,
                                                         annotations,genomeSequence, stranded = stranded, verbose = verbose)
        # TODO: [OTHER] runName = "TODO" is a placeholder; replace with a meaningful sample/run identifier
        se <- isore.constructReadClasses(readGrgList,
                                              unlisted_junctions, uniqueJunctions, runName = "TODO",
                                              annotations, stranded, verbose)

    }
    metadata(se)$warnings <- warnings
    if(trackReads){
        metadata(se)$readNames <- names(readGrgList)
        metadata(se)$readId <- mcols(readGrgList)$id
    }
    rm(readGrgList)
    refSeqLevels <- seqlevels(genomeSequence)
    GenomeInfoDb::seqlevels(se) <- refSeqLevels
    # create SE object with reconstructed readClasses
    se <- scoreReadClasses(se, genomeSequence, annotations, 
                             defaultModels = defaultModels,
                             fit = fitReadClassModel,
                             returnModel = returnModel,
                             min.readCount = min.readCount,
                             min.exonOverlap = min.exonOverlap,
                             fusionMode = fusionMode,
                             verbose = verbose)
    return(se)
}


# --- lowMemoryConstructReadClasses ---
# Module: Module 2 — Read processing (per sample) | bambu-processReads.R
# Called by: bambu-processReads.R
# Call count: 2 calls, 1 files
#' Low memory mode for construct read classes (processByChromosome)
#' @noRd
lowMemoryConstructReadClasses <- function(readGrgList, genomeSequence, 
                                          annotations, stranded, verbose,bam.file, fusionMode = FALSE){
    if(fusionMode){
        readGrgList <- list(readGrgList)
        names(readGrgList) <- c("fusion")
    } else{
        readGrgList <- split(readGrgList, getChrFromGrList(readGrgList))
    }
    se <- lapply(names(readGrgList),FUN = function(i){
        if(length(readGrgList[[i]]) == 0) return(NULL)
        # create error and strand corrected junction tables
        unlisted_junctions <- unlistIntrons(readGrgList[[i]], use.ids = TRUE)
        uniqueJunctions <- isore.constructJunctionTables(unlisted_junctions,
                                                         annotations,genomeSequence, stranded = stranded, verbose = verbose)
        # TODO: [OTHER] runName = "TODO" is a placeholder; replace with a meaningful sample/run identifier (e.g. chromosome name i)
        se.temp <- isore.constructReadClasses(readGrgList[[i]],
                                              unlisted_junctions, uniqueJunctions, runName = "TODO",
                                              annotations, stranded, verbose)
        return(se.temp)
    })
    se <- se[!sapply(se, FUN = is.null)]
    se <- do.call("rbind",se)
    rownames(se) <- paste("rc", seq_len(nrow(se)), sep = ".")
    return(se)
}

# --- seqlevelCheckReadsAnnotation ---
# Module: Module 2 — Read processing (per sample) | bambu-processReads.R
# Called by: bambu-processReads.R
# Call count: 2 calls, 1 files
#' Check seqlevels for reads and annotations
#' @importFrom GenomeInfoDb seqlevels
#' @noRd
seqlevelCheckReadsAnnotation <- function(reads, annotations){ #TODO (JG) [validate-input] should this be done by reading annotations and bam file for chromosome style match? instead of here? downstream of checkinput section we can assume all is correc then?
    warnings <- c()
    if (length(intersect(seqlevels(reads),
                         seqlevels(annotations))) == 0)
        warnings <- c(warnings, paste0("no annotations with matching seqlevel styles, ",
        "all missing chromosomes will use de-novo annotations"))
    if (!all(seqlevels(reads) %in% 
             seqlevels(annotations))) 
        warnings <- c(warnings, paste0("not all chromosomes present in reference annotations, ",
            "annotations might be incomplete. Please compare objects ",
            "on the same reference"))
    return(warnings)
}


#' Split read class files
#' @importFrom dplyr Matrix
#' @noRd
splitReadClassFiles = function(readClassFile){  #TODO (JG) [bambu-modules] this is only used in assignDist, move there? what is this function doing
    distTable <- metadata(metadata(readClassFile)$readClassDist)$distTable  
    eqClasses <- distTable %>% group_by(eqClassById) %>% 
        distinct(eqClassById, readCount,GENEID, totalWidth, firstExonWidth, .keep_all = TRUE)
    eqClasses$sampleIDs <- rowData(readClassFile)$sampleIDs[match(eqClasses$readClassId, rownames(readClassFile))]
    eqClasses <- eqClasses %>% summarise(nobs = sum(readCount),
                                                sampleIDs = list(unlist(sampleIDs)))
    counts.table <- tableFunction(eqClasses$sampleIDs)
    counts <- sparseMatrix(
        i = rep(seq_along(counts.table), lengths(counts.table)),
        j = as.numeric(names(unlist(counts.table))),
        x = unlist(counts.table),
        dims = c(nrow(eqClasses), nrow(metadata(readClassFile)$sampleData)))
    #incompatible counts
    distTable <- metadata(metadata(readClassFile)$readClassDist)$distTable.incompatible
    if(nrow(distTable)==0) {
        counts.incompatible <- sparseMatrix(i= 1, j = 1, x = 0,
        dims = c(1, length(metadata(readClassFile)$sampleData$id)))
        rownames(counts.incompatible) <- "TODO"
    } else{
        distTable$sampleIDs <- rowData(readClassFile)$sampleIDs[match(distTable$readClassId, rownames(readClassFile))]
        distTable <- distTable %>% group_by(GENEID.i) %>% summarise(counts = sum(readCount),
                    sampleIDs = list(unlist(sampleIDs)))
        counts.table <- lapply(distTable$sampleIDs, FUN = function(x){table(x)})
        counts.incompatible <- sparseMatrix(
            i = rep(seq_along(counts.table), lengths(counts.table)),
            j = as.numeric(names(unlist(counts.table))),
            x = unlist(counts.table),
            dims = c(nrow(distTable), length(metadata(readClassFile)$sampleData$id)))
        colnames(counts.incompatible) <- metadata(readClassFile)$sampleData$id
        rownames(counts.incompatible) <- distTable$GENEID.i 
    }
    colnames(counts) <- metadata(readClassFile)$sampleData$id
    metadata(readClassFile)$eqClassById <- eqClasses$eqClassById
    #rownames(counts) = eqClasses$eqClassById
    metadata(readClassFile)$countMatrix <- counts
    metadata(readClassFile)$incompatibleCountMatrix <- counts.incompatible  
    return(readClassFile)
}


# --- splitReadClassFilesByRC ---
# Module: Module 2 — Read processing (per sample) | bambu-processReads.R
# Called by: bambu-extendAnnotations-utilityExtend.R
# Call count: 1 calls, 1 files
#' Split read class files by RC
#' @importFrom Matrix
#' @noRd
splitReadClassFilesByRC <- function(readClassFile){ #TODO (JG) [bambu-modules] this is only used in bambu, part of clustering. Should not be here
    counts.table <- tableFunction(rowData(readClassFile)$sampleIDs)
    counts <- sparseMatrix(
        i = rep(seq_along(counts.table), lengths(counts.table)),
        j = as.numeric(names(unlist(counts.table))),
        x = unlist(counts.table),
        dims = c(nrow(readClassFile), length(metadata(readClassFile)$samples)))
    return(counts)
}

# --- tableFunction ---
# Module: Module 2 — Read processing (per sample) | bambu-processReads.R
# Called by: bambu-processReads.R
# Call count: 2 calls, 1 files
#' table sample IDs list column
#' @noRd
tableFunction <- function(xList){ #TODO (JG) [bambu-modules] this is only used as part of clustering, should move with the above function
    return(lapply(xList, function(x) table(x)))
}
