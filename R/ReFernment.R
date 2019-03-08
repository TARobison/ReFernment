#' @title Annotate RNA editing in for plastid genome submission
#'
#' @description Provides RNA editing annotations to satisfy GenBank submission requirements. 
#' Takes as input the paths to the folders containing required file inputs and
#' the desired output path and the names of 
#' the files that the user wishes to annotate. Each filetype needn't be in 
#' separate folders, but it is recommended. This is especially true for output 
#' file path, because if it isn't separate from input gb files, the original 
#' files will be overwritten! Be sure that the file names for a respective 
#' genome are the same across file types! 
#' 
#' @param gbFolderPath this is the path to the folder containing the gb file(s)
#' that you would like annotated. 
#' @param gffFolderPath this is the path to the folder containing the gff 
#' file(s) that you would like annotated. 
#' @param fastaFolderPath this is the path to the folder containing the 
#' fasta file(s) corresponding to the gff/gb files you would like annotated.
#' @param outFolderPath this is the path to the folder where the new genome
#' annotations will be produced. In addition to the preferred output file (gb or gff), the new 
#' annotations will produced as a feature table and a protein fasta file for BankIt submission. 
#' It is important to note that you
#' should have the output file separate from your input files, as this prevents
#' the script from overwriting the original files!
#' @param genomes this should be a list containing the filenames of all the 
#' genomes that you would like to annotate. Note that you should not include
#' the file endings of these files. 
#'
#' @return an annotated gb or gff file file. In addition, a five column feature table will be produced, which can be used for BankIt submissions. This file will have the required misc_feature
#' annotations as well as a translation corrected for RNA editing.   
#'
#' @author Tanner Robison
#' @examples
#' 
#' genomes <- c("Asplenium_pek", "Woodwardia_uni")
#' gbFolder <- "/home/travis/build/TARobison/ReFernment/examples/GB/"
#' gffFolder <- "/home/travis/build/TARobison/ReFernment/examples/GFF/"
#' outputFolder <- "/home/travis/build/TARobison/ReFernment/examples/"
#' ReFernment(gffFolder, fastaFolderPath=NULL, gbFolderPath=NULL, outputFolder, genomes)
#' 
#' @importFrom ape read.gff
#' @importFrom Biostrings DNAString DNAStringSet readDNAStringSet reverseComplement translate reverseComplement getGeneticCode
#' @importFrom stringr str_extract str_extract_all str_to_lower
#' @importFrom utils write.table read.csv write.csv
#' @importFrom stats na.omit
#' @importFrom dplyr arrange
#' 
#' @export


ReFernment <- function(gffFolderPath, fastaFolderPath=NULL, gbFolderPath=NULL, outFolderPath, genomes){
    
    editGB_File <- function(genomeName){
        #this function acts as a wrapper for many of the other functions contained within ReFernment. It adds RNA editing annotations as well as conceptual translations for the gb file. 
        #the protein fasta file is also made here. 
        proteinFasta <- matrix(NA, 1,4)
        for(i in 1:nrow(annotations)){
            if(annotations$type[i] == "CDS"){
                intervalProblem <- FALSE
                gene <- getGeneName(annotations[i,])
                exonIndex <- getNumberOfExons(annotations, i)
                exonCount <- length(exonIndex)
                #cat(gene)
                if(identical(annotations$attributes[i], annotations$attributes[i-1]) && identical(annotations$strand[i], annotations$strand[i-1])){next}
                else if(identical(annotations$attributes[i], annotations$attributes[i+1]) && identical(annotations$strand[i], annotations$strand[i+1])){
                    for(j in 0:exonCount-1){
                        if(annotations$start[i+j] > nchar(dnachar) | annotations$end[i+j] > nchar(dnachar)){
                            intervalProblem <- TRUE
                            break
                        }
                    }
                    if(isTRUE(intervalProblem)){
                        cat("The intervals for", gene, "'s annotation are incorrect. Part or all of the annotation extend beyond the sequence length. Please check.\n")
                        next
                    }
                    exons <- checkExons(annotations,i, dnachar,gb)  
                    gb <- exons[[1]]
                    protein <- exons[[2]]
                    annotations <- exons[[3]]
                    proteinProduct <- getGeneProduct(annotations[i,])
                    proteinrow <- cbind(gene, protein, annotations$start[i], proteinProduct)
                    names(proteinrow) <- names(proteinFasta)
                    proteinFasta <- rbind(proteinFasta, proteinrow)
                }else{
                    
                    CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
                    if(nchar(CDS) %% 3 != 0){
                        outCDS <- checkFrame(annotations, i, CDS,gb)
                        gb <- outCDS[[1]]
                        annotations <- outCDS[[2]]
                        CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
                    }
                     
                    outStart <- checkStartCodon(CDS, annotations,i,gb)
                    gb <-outStart[[1]]
                    annotations <-outStart[[2]]
                    CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])

                    outStop <- checkStopCodon(CDS, annotations,i, gb)
                    gb <- outStop[[1]]
                    annotations <- outStop[[2]]
                    CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])

                    outInternal <- checkForInternalStops(CDS, annotations, i, gb)
                    gb <- outInternal[[1]]
                    annotations <- outInternal[[2]]
                    
                    translation <- correctTranslation(annotations, i, dnachar, CDS)
                    protein <- translation[1]
                    proteinProduct <- getGeneProduct(annotations[i,])
                    proteinrow <- cbind(gene, protein, annotations$start[i], proteinProduct)
                    names(proteinrow) <- names(proteinFasta)
                    RNAediting <- translation[2]
                    gb <- annotateTranslation(protein, RNAediting, annotations, gb,i)
                    proteinFasta <- rbind(proteinFasta, proteinrow)
                  #  cat(gene)
                    }
            }
        }
        writeLines(gb, paste(output, ".gb", sep=''))
        write.csv(proteinFasta, paste(output, ".csv", sep=''), row.names=FALSE)
    }
    
    edit_GFF_file <-function(genomeName){
        proteinFasta <- matrix(NA, 1,4)
        for(i in 1:nrow(annotations)){
            intervalProblem <- FALSE
            if(annotations$type[i] == "CDS"){
                gene <- getGeneName(annotations[i,])
                exonIndex <- getNumberOfExons(annotations, i)
                exonCount <- length(exonIndex)
                if(identical(annotations$attributes[i], annotations$attributes[i-1]) && identical(annotations$strand[i], annotations$strand[i-1])){next}
                else if(identical(annotations$attributes[i], annotations$attributes[i+1]) && identical(annotations$strand[i], annotations$strand[i+1])){
                    for(j in 0:exonCount){
                        if(annotations$start[i+j] > nchar(dnachar) | annotations$end[i+j] > nchar(dnachar)){
                            intervalProblem <- TRUE
                            break
                        }
                    }
                    if(isTRUE(intervalProblem)){
                        cat("The intervals for", gene, "'s annotation are incorrect. Part or all of the annotation extend beyond the sequence length. Please check.\n")
                        next
                    }
                    exons <- checkExons(annotations,i, dnachar,gb)  
                    gb <- exons[[1]]
                    protein <- exons[[2]]
                    annotations <- exons[[3]]
                    #yell at people if they don't have a product in their gff
                    proteinProduct <- getGeneProduct(annotations[i,])
                    proteinrow <- cbind(gene, protein, annotations$start[i], proteinProduct)
                    names(proteinrow) <- names(proteinFasta)
                    proteinFasta <- rbind(proteinFasta, proteinrow)
                }else{
                    
                    CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
                    if(nchar(CDS) %% 3 != 0){
                        outCDS <- checkFrame(annotations, i, CDS,gb)
                        gb <- outCDS[[1]]
                        annotations <- outCDS[[2]]
                        CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
                    }

                    outStart <- checkStartCodon(CDS, annotations,i,gb)
                    gb <-outStart[[1]]
                    annotations <-outStart[[2]]
                    CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])

                    outStop <- checkStopCodon(CDS, annotations,i, gb)
                    gb <- outStop[[1]]
                    annotations <- outStop[[2]]
                    CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])

                    outInternal <- checkForInternalStops(CDS, annotations, i, gb)
                    gb <- outInternal[[1]]
                    annotations <- outInternal[[2]]
                    
            
                    translation <- correctTranslation(annotations, i, dnachar, CDS)
                    protein <- translation[1]
                    proteinProduct <- getGeneProduct(annotations[i,])
                    proteinrow <- cbind(gene, protein, annotations$start[i], proteinProduct)
                    names(proteinrow) <- names(proteinFasta)
                    RNAediting <- translation[2]
                    gb <- annotateTranslation(protein, RNAediting, annotations, gb,i)
                    proteinFasta <- rbind(proteinFasta, proteinrow)
                    }
            }
        }
    
        write.table(annotations, paste(output, ".gff", sep=''), quote=FALSE, sep='\t', na='.', row.names=FALSE, col.names=FALSE)
        write.csv(proteinFasta, paste(output, ".csv", sep=''), row.names=FALSE)
        return(annotations)
    }   
    
    checkFrame <- function(annotations, i, CDS, gb){
        gene <- getGeneName(annotations[i[1],])
        numExons <- length(i)
        if(numExons== 1){
            CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
        }   
        if(annotations$strand[i[numExons]] == '-'){
            if(!is.null(gb)){
                original <- annotations$start[i[numExons]]
                while((nchar(CDS) %% 3) != 0){
                    annotations$start[i[numExons]] <- annotations$start[i[numExons]] - 1
                    CDS <- paste(CDS, "N", sep='')
                }
                gbString <- paste(original, '(\\.\\.)', sep='')
                newStop <- paste(annotations$start[i[numExons]],'..', sep='')
                gb <- gsub(gbString, newStop,gb)
                out <- list(gb, annotations)
                return(out)
            }else{
                while((nchar(CDS) %% 3) != 0){
                    annotations$start[i[numExons]] <- annotations$start[i[numExons]] - 1
                    CDS <- paste(CDS, "N", sep='')
                }
                return(annotations)
            }
        }
        if(annotations$strand[i[numExons]] == '+'){
            if(!is.null(gb)){
                original <- annotations$end[i[numExons]]
                while((nchar(CDS) %% 3) != 0){
                    annotations$end[i[numExons]] <- annotations$end[i[numExons]] +1
                    CDS <- paste(CDS, "N", sep='')
                }
                gbString <- paste('(\\.\\.)', original, sep='')
                newStop <- paste('..',annotations$end[i[numExons]], sep='')
                gb <- gsub(gbString, newStop,gb)
                out <- list(gb, annotations)
                return(out)
            }
        }else{
            while((nchar(CDS) %% 3) != 0){
                annotations$end[i[numExons]] <- annotations$end[i[numExons]] +1
                CDS <- paste(CDS, "N", sep='')
            }
            return(annotations)
        }
    }
    
    concatenateExons <- function(numExons, exonIndex, annotations, dnachar, i){
        #concatenates exons so we have a single coding sequence.
        firstExon <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
        for(j in 2:numExons){
                nextExon <- getCodingSequence(annotations[exonIndex[j],], dnachar, annotations$start[exonIndex[j]], annotations$end[exonIndex[j]])
                firstExon <- paste(firstExon, nextExon, sep = '')
        }
        return(firstExon)
    }

    checkExons <- function(annotations,i, dnachar, gb){
        #a wrapper function for genes with multiple intervals.
        exonIndex <- getNumberOfExons(annotations, i)
        numExons <- length(exonIndex)
        CDS <- concatenateExons(numExons, exonIndex, annotations, dnachar, i)
    
        if(nchar(CDS) %% 3 != 0){
            if(is.null(gb)){
                annotations <- checkFrame(annotations, i, CDS, gb)
            }else{
                outCDS <- checkFrame(annotations, i, CDS,gb)
                gb <- outCDS[[1]]
                annotations <- outCDS[[2]]
                CDS <- concatenateExons(numExons, exonIndex, annotations, dnachar, i)
            }
        }
        

        ## THIS IS A POTENTIAL DANGER ZONE, but I think it will be fine if I am careful with the output of the functions below
        outStart <- checkStartCodon(CDS, annotations, exonIndex[1], gb)
        gb <- outStart[[1]]
        annotations <- outStart[[2]]
        CDS <- concatenateExons(numExons, exonIndex, annotations, dnachar, i)

        finalExon <- getCodingSequence(annotations[exonIndex[numExons],], dnachar, annotations$start[exonIndex[numExons]], annotations$end[exonIndex[numExons]])
        outStop <- checkStopCodon(finalExon, annotations,exonIndex, gb) ###removed exonIndex[numExons]
        gb <- outStop[[1]]
        annotations <- outStop[[2]]
        CDS <- concatenateExons(numExons, exonIndex, annotations, dnachar, i)

        outInternal <- checkForInternalStops(CDS, annotations,exonIndex,gb)
        gb <- outInternal[[1]]
        annotations <- outInternal[[2]]


        translation <- correctTranslation(annotations, i, dnachar,CDS)
        protein <- as.character(translation[1])
        RNAediting <- translation[2]
        gb <- annotateTranslation(protein, RNAediting, annotations,  gb,i)
        return(list(gb, protein, annotations))
    }

    correctTranslation <- function(annotations, i, dnachar, CDS){
        #corrects coding sequences that would be restored by RNA editing. 

        RNAediting = FALSE
        
        exonIndex <- getNumberOfExons(annotations, i)

        if(length(exonIndex) > 1){
            codonList <- codonGroup(CDS)
        }else{
        codonList <- codonGroup(getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i]))
        }
        if(codonList[1] == possibleEditedStarts){
            codonList[1] <- CtoU_RNAediting(codonList[1])
            RNAediting = TRUE
        }
        if(any(codonList[length(codonList)] == possibleEditedStops)){
            codonList[length(codonList)] <- UtoC_RNAediting(codonList[length(codonList)])
            codonList[length(codonList)] <- CtoU_RNAediting(codonList[length(codonList)])
            RNAediting = TRUE
        }
        for(j in 2:length(codonList)-1){
            if(any(codonList[j] == conventionalStops)){
                codonList[j] <-  UtoC_RNAediting(codonList[j])
                RNAediting = TRUE
            }
        }

        correctedSeq <- paste(codonList, collapse='')
       
        correctedTranslation <- as.character(translate(DNAStringSet(correctedSeq), getGeneticCode("11"),if.fuzzy.codon="solve"))
        
        out <- list(correctedTranslation, RNAediting)

        return(out)
    }

    annotateTranslation <- function(correctedTranslation, RNAediting, annotations, gb, i){
        #edits GB file to add RNA editing qualifiers as well as conceptual translations.
        
        exonIndex <- getNumberOfExons(annotations, i)
        
        if(is.null(gb)){
            return(gb)
        }else{
            if(isTRUE(as.logical(RNAediting))){
                gbInsert <- paste('/translation="',correctedTranslation, '"\n                     /exception="RNA editing"',sep='')
            }else{
                gbInsert <- paste('/translation="',correctedTranslation, '"',sep='')
            }
        
            if(length(exonIndex) > 1){
                if(length(exonIndex) >2){
                    if(annotations[i[1], "strand"] == '-'){
                        gbstring<- paste('(\\s{21}complement\\(.*',annotations[exonIndex[length(exonIndex)],"start"],'.*)',sep='')  #annotations[exonIndex[length(exonIndex)],"start"]
                    }
                    if(annotations[i[1],"strand"] == '+'){
                        # gbstring<-paste('(CDS\\s*join\\(',annotations[i,"start"],'.*)',sep='')
                        gbstring<- paste('(\\s{21}.*',annotations[exonIndex[length(exonIndex)],"start"],'.*)',sep='')
                    }
                }else{
                    if(annotations[i[1], "strand"] == '-'){
                        gbstring<- paste('(CDS\\s*join\\(complement\\(',annotations[exonIndex[1],"start"],'.*)',sep='')  #annotations[exonIndex[length(exonIndex)],"start"]
                        additionalLine <- grep(paste('(\\s{21}complement\\(', annotations$start[exonIndex[2]], ')', sep =''), gb)
                        if(length(additionalLine) != 0){
                        gbstring <- paste('(\\s{21}complement\\(', annotations$start[exonIndex[2]], '.*)', sep ='')
                        }   
                    }
                    if(annotations[i[1],"strand"] == '+'){
                        gbstring<-paste('(CDS\\s*join\\(',annotations[i,"start"],'.*)',sep='')
                    
                        # additionalLine <- grep(paste('(\\s{21}complement\\(', annotations$start[exonIndex[2]], ')', sep =''), gb)
                        # if(length(additionalLine) != 0){
                        # gbstring <- paste('(\\s{21}complement\\(', annotations$start[exonIndex[2]], '.*)', sep ='')
                    }
                }
            }else{
                if(annotations[i[1], "strand"] == '-'){
                    gbstring<- paste('(CDS\\s*complement\\(',annotations[i,"start"],'\\.\\.',annotations[i,"end"],'\\))',sep='')
                    additionalLine <- grep(paste('(\\s{21}complement\\(', annotations$start[exonIndex[2]], ')', sep =''), gb)

                    if(length(additionalLine) != 0){
                        gbstring <- paste('(\\s{21}complement\\(', annotations$start[exonIndex[2]], '.*)', sep ='')
                    }
                }
                if(annotations[i[1],"strand"] == '+'){
                    gbstring<-paste('(CDS\\s*',annotations[i,"start"],'\\.\\.',annotations[i,"end"],')',sep='')
                }
            }
            newString <-  paste("\\1\n\                     ",gbInsert, sep='')
            gb <- sub(gbstring[1], newString[1], gb)
            
            return(gb)
        }
    }

    checkStartCodon <- function(CDS, annotations,i,gb){
        #checks the first residue to see if it is a valid start codon. If not, checks to see whether RNA editing, or changing the gene boundaries restores start codon
        originalGFF <- annotations
        if(any(substr(CDS, 1, 3) == conventionalStarts)){
            out <- list(gb, annotations)
            return(out)
        }
        possibleEditedStarts <- c('ACG')
        if(any(substr(CDS, 1, 3) == possibleEditedStarts)){
            #CHECK THE ADD FEATURE FUNCTION. MAKE A GB FILE LOGICAL
            if(annotations[i[1], "strand"] == '+'){
                out <- addFeature(annotations, i, annotations[i, "start"]+1, "CtoU", gb, "gene")
            }
            if(annotations[i[1], "strand"] == '-'){
                
                out <- addFeature(annotations, i, annotations[i, "end"]-1, "CtoU",gb,"gene")
            }
            return(out)
        }else{
            shortenOut <- shortenUpstream(annotations, dnachar, gb, CDS,i)
            gb <- shortenOut[[1]]
            annotations <- shortenOut[[2]]
            if(identical(originalGFF, annotations)){
                #####
                out <- extendUpstream(annotations, dnachar, gb, i)
                return(out)
            }
            out <- list(gb, annotations)
            return(out)
        }   
    }

    checkStopCodon <- function(CDS, annotations,i, gb){
        #checks the final residue to see if it is a valid stop codon. If not, checks to see whether RNA editing, or changing the gene boundaries restores stop codon.
        gene <- getGeneName(annotations[i[1],])
        exonIndex <- getNumberOfExons(annotations, i[1])
        numExons <- length(exonIndex)
        codonList <- codonGroup(CDS)
        originalGFF <- annotations
        if(any(substr(CDS, nchar(CDS)-2,nchar(CDS)) == (conventionalStops))){
            out <- list(gb, annotations)
            return(out) 
        } else if(any(substr(CDS, nchar(CDS)-2,nchar(CDS)) == possibleEditedStops)){
            if(annotations$strand[i[numExons]] == '+'){
                out <- addFeature(annotations, i[1], annotations[i[numExons], "end"]-2, "CtoU",gb,"gene") #### changed these
            }
            if(annotations$strand[i[numExons]] == '-'){
                out <- addFeature(annotations, i[1],annotations[i[numExons], "start"]+2, "CtoU", gb, "gene") #######
            }
            return(out) 
        } else {
            shortenOut <- shortenDownstream(annotations, dnachar, gb, CDS, i)
            gb <- shortenOut[[1]]
            annotations <- shortenOut[[2]]
            if(identical(originalGFF, annotations)){
                #####
                out <- extendDownstream(annotations, dnachar, gb, i)
                return(out)
            }
            out <- list(gb, annotations)
            return(out)
        }
    }

    extendDownstream <- function(annotations, dnachar, gb, i){
        #if the stop codon is not valid, this function may extend the downstream end of gene to check if valid stop lies nearby
        if(annotations$strand[i[1]] == '-'){
            CDS <- getCodingSequence(annotations[i,], dnachar, (annotations$start[i]-30), annotations$end[i])
        }
        if(annotations$strand[i[1]] == '+'){
            CDS <- getCodingSequence(annotations, dnachar, annotations$start[i], annotations$end[i] +30)
        }
        codonList <- codonGroup(CDS)
        exonIndex <- getNumberOfExons(annotations, i)
        
        numExons <- length(exonIndex)
        for(j in (length(codonList)-10):((length(codonList)))){
            if(any(codonList[j] == conventionalStops)){
                if(annotations$strand[i[numExons]] == '-'){
                    dist <- j - (length(codonList)-10)
                    newStop <- annotations$start[i[numExons]] - dist*3
                    if(is.null(gb)){
                        annotations$start[i[numExons]] <- newStop
                    }else{
                        gbstring <- paste(annotations$start[i[numExons]],'\\.\\.',sep='')
                        gb <- sub(gbstring, paste(newStop,'\\.\\.',sep=''), gb, perl = TRUE)
                        annotations$start[i[numExons]] <- newStop
                    }
                    out <- list(gb, annotations)
                    return(out)
                }
                if(annotations$strand[i[numExons]] == '+'){
                    dist <- j -(length(codonList) -10)
                    newStop <- annotations$end[i[numExons]] + dist*3
                    if(is.null(gb)){
                        annotations$end[i[numExons]] <- newStop
                    }else{
                        gbstring <- paste('\\.\\.',annotations$end[i[numExons]],sep='')
                        gb <- sub(gbstring, paste('\\.\\.', newStop,sep='') , gb, perl = TRUE)
                        annotations$end[i[numExons]] <- newStop
                    }
                    out <- list(gb, annotations)
                    return(out)
                }  
            }
        }
        out <- list(gb, annotations)
        return(out)   
    }
 
    extendUpstream <- function(annotations,dnachar, gb, i){
        #if the start codon is not valid, this function may extend the upstream end of gene to check if valid start codon lies nearby
        if(annotations$strand[i[1]] == '-'){
            CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i]+18)
        }
        if(annotations$strand[i[1]] == '+'){
            CDS <- getCodingSequence(annotations, dnachar, annotations$start[i] -18, annotations$end[i])
        }
        codonList <- codonGroup(CDS)
        count <- 0
        for(j in 6:1){
            count <- count+1
            if(any(codonList[j] == conventionalStarts)){
                if(annotations$strand[i[1]] == '-'){
                    newStart <- annotations$end[i[1]] + count*3
                    if(is.null(gb)){
                        annotations$end[i[1]] <- newStart
                    }else{
                        gbstring <- paste('\\.\\.',annotations$end[i[1]], '\\)',sep='')
                        gb <- sub(gbstring, paste('\\.\\.', newStart[1], ')', sep=''), gb, perl = TRUE)
                        annotations$end[i[1]] <- newStart
                    }
                    out <- list(gb, annotations)
                    return(out)
                }
                if(annotations$strand[i[1]] == '+'){
                    newStart <- annotations$start[i[1]] - count*3
                    if(is.null(gb)){
                        annotations$start[i[1]] <- newStart
                    }else{
                        gbstring <- paste(annotations$start[i[1]],'\\.\\.',sep='')
                        gb <- sub(gbstring, paste(newStart,'\\.\\.',sep='') , gb, perl = TRUE)
                        annotations$start[i[1]] <- newStart
                    }
                    out <- list(gb, annotations)
                    return(out)
                }  
            }
        }
        out <- list(gb, annotations)
        return(out)   
    }

    checkForInternalStops <- function(CDS, annotations, i, gb){
        #moves through each residue of a gene. If a stop codon is detected, function checks whether gene could be restored through RNA editing. 
        gene <- getGeneName(annotations[i,])
        codonList <-codonGroup(CDS)
        stops <- 0
        totalLength <- 0
        exonCount <- 1
        out <- list(gb, annotations)
            for(j in 1:(length(codonList)-1)){            
                if(annotations$strand[i[1]] == '-'){
                    if(annotations$end[i[exonCount]] - annotations$start[i[exonCount]] +totalLength < j*3 && exonCount < length(i)){ 
                        #this is the most ridiculous way to fix this problem. Could barely read it today (maybe I needed more coffee though)
                        #fix in the future
                        totalLength <- annotations$end[i[exonCount]] - annotations$start[i[exonCount]] + totalLength +1
                        exonCount <- exonCount +1
                    }
                }
                if(annotations$strand[i[1]] == '+'){
                    if(annotations$end[i[exonCount]] - annotations$start[i[exonCount]] +totalLength < j*3 && exonCount < length(i)){
                        totalLength <- annotations$end[i[exonCount]] - annotations$start[i[exonCount]] + totalLength +1
                        exonCount <- exonCount +1
                    }
                }
                if(any(codonList[j] == conventionalStops)){
                    if(annotations$strand[i[exonCount]] == '-'){
                        termStart <- annotations$end[i[exonCount]] - ((j-1)*3 - totalLength)
                        #changed for new output
                        out <-addFeature(annotations, i[1], termStart,"UtoC",gb,gene)
                        annotations <- out[[2]]
                        gb <- out[[1]]
                        stops <- stops +1
                    }else{
                        termStart <- annotations$start[i[exonCount]] + ((j-1)*3 - totalLength)
                        #changed for new output
                        out <-addFeature(annotations, i[1], termStart,"UtoC",gb,gene)
                        annotations <- out[[2]]
                        gb <- out[[1]]
                        stops <- stops +1
                    }

                } 
            }
        if(stops > 5){
            gene <- getGeneName(annotations[i,])
            cat("There are a high number of edited Stops (", stops,") in", gene, "manually check to make sure frame is correct\n")
            #browser()
        }
        return(out)
    }

    addFeature <- function(annotations, i, site, editType, gb, gene){
        #adds RNA editing qualifiers to gb file. 
        exonIndex <- getNumberOfExons(annotations,i)
        if(is.null(gb)){
            newFeature <- annotations[nrow(annotations),]
            if(editType == 'CtoU'){
                newFeature$attributes <- "Name=C to U RNA editing;note=putative C to U RNA editing"
            }
            if(editType == 'UtoC'){
                newFeature$attributes <- "Name=U to C RNA editing;note=putative U to C RNA editing"
            }
            newFeature$type <- "misc_feature"
            newFeature$start <- site
            newFeature$end <- site
            newFeature$strand <- NA
            annotations <- rbind(annotations, newFeature)
            rownames(annotations) <- 1:nrow(annotations)
        }else{
            if(editType == 'CtoU'){
                note <- '/note="putative C to U RNA editing"'
            }
            if(editType == 'UtoC'){
                note <- '/note="putative U to C RNA editing"'
            }
            newFeature <- paste('     misc_feature    ', site, '\n\                     ', note,'\n', sep ='')
            if(length(exonIndex) > 1){
                if(annotations[i[1], "strand"] == '-'){
                    gbstring<- paste('(     CDS\\s*join\\(complement\\(',annotations[i,"start"],'.*)',sep='')
                }
                if(annotations[i[1],"strand"] == '+'){
                    gbstring<- paste('(     CDS\\s*join\\(',annotations[i,"start"][1],'.*)',sep='')
                }
            }else{
                if(annotations[i[1], "strand"] == '-'){
                    gbstring<- paste('(     CDS\\s*complement\\(',annotations[i,"start"],'\\.\\.',annotations[i,"end"],'\\))',sep='')
                }
                if(annotations[i[1],"strand"] == '+'){
                    gbstring<-paste('(     CDS\\s*',annotations[i,"start"],'\\.\\.',annotations[i,"end"],')',sep='')
                }
            }
            newString <-  paste(newFeature, "\\1",sep='')
            gb <- sub(gbstring[1], newString[1], gb)
            
        }
        out  <- list(gb, annotations)
        return(out)
    }

    shortenDownstream <- function(annotations, dnachar, gb, CDS,i){
        #if the stop codon is not valid, this function may shorten the downstream end of gene to check if valid stop lies nearby        
        codonList <-codonGroup(CDS)
        
        exonIndex <- getNumberOfExons(annotations, i)
        numExons <- length(exonIndex)
        for(j in (length(codonList)-7):((length(codonList)-1))){
            if(any(codonList[j] == conventionalStops)){
                if(annotations$strand[i[numExons]] == '-'){
                    dist <- length(codonList) - j
                    newStop <- annotations$start[i[numExons]] + dist*3
                    if(!is.null(gb)){
                        gbstring <- paste(annotations$start[i[numExons]],'\\.\\.',sep='')
                        gb <- sub(gbstring, paste(newStop,'\\.\\.',sep=''), gb, perl = TRUE)
                    }
                    annotations$start[i[numExons]] <- newStop
                    out <- list(gb, annotations)
                    return(out)
                }
                if(annotations$strand[i[numExons]] == '+'){
                    dist <- length(codonList) - j
                    newStop <- annotations$end[i[numExons]] - dist*3
                    if(!is.null(gb)){
                        gbstring <- paste('\\.\\.',annotations$end[i[numExons]],sep='')
                        gb <- sub(gbstring, paste('\\.\\.', newStop,sep='') , gb, perl = TRUE)
                    }
                    annotations$end[i[numExons]] <- newStop
                    out <- list(gb, annotations)
                    return(out)
                }  
            }
        }
        out <- list(gb, annotations)
        return(out)
    }

    shortenUpstream <- function(annotations, dnachar, gb, CDS,i){
        #if the start codon is not valid, this function may shorten the upstream end of gene to check if valid start codon lies nearby
        codonList <-codonGroup(CDS)
        for(j in 1:5){
            if(any(codonList[j] == conventionalStarts)){
                if(annotations[i[1],"strand"] == '-'){
                    newStart <- annotations[i,"end"] - (j-1)*3
                    if(is.null(gb)){
                        gbstring <- paste('\\.\\.',annotations[i,"end"],sep='')
                        gb <- sub(gbstring, paste('\\.\\.', newStart[1],sep=''), gb, perl = TRUE)
                    }
                    annotations[i,"end"] <- newStart
                    out <- list(gb, annotations)
                    return(out)
                }
                if(annotations[i[1],"strand"] == '+'){
                    newStart <- annotations[i,"start"] + (j-1)*3
                    if(is.null(gb)){
                        gbstring <- paste(annotations[i,"start"], '\\.\\.',sep='')
                        gb <- sub(gbstring, paste(newStart, '\\.\\.',sep='') , gb, perl = TRUE)
                    }
                    annotations[i,"start"] <- newStart
                    out <- list(gb, annotations)
                    return(out)
                }    
            }
        }
        out <- list(gb, annotations)
        return(out)   
    }

    GFF2FeatureTable <- function(annotations, outPath, name){
        illegalQualifiers <- c("NCBI Feature Key", "NCBI Join Type", "ID", "codon_start", "Name", "mol_type", "db_xref", "translation") 
        illegalFeatures <- c("region", "source")   
        featureTable <- matrix('', 1, 5)
        featureTable[1,1] <- paste(">Feature", name)
       
        for(i in 1:nrow(annotations)){
            if(any(str_to_lower(annotations$type[i]) == str_to_lower(illegalFeatures))){next}
            gene <- getGeneName(annotations[i,])
            newRow <- matrix('', 1, 5)
            if(identical(gene, getGeneName(annotations[i-1,]))
                    && identical(annotations$strand[i], annotations$strand[i-1])
                    && !is.na(annotations$strand[i])){
                newRow[1,3] <- ''
            }else{
                newRow[1,3] <- as.character(annotations$type[i])
            }
            if(is.na(annotations$strand[i])){
                newRow[1,1] <- annotations$start[i]
                newRow[1,2] <- annotations$end[i]
                featureTable <- rbind(featureTable, newRow)
            }else if(annotations$strand[i] == '-'){
                newRow[1,1] <- annotations$end[i]
                newRow[1,2] <- annotations$start[i]
                featureTable <- rbind(featureTable, newRow)
                if(identical(gene, getGeneName(annotations[i+1,]))
                    && identical(annotations$strand[i], annotations$strand[i+1])){next}
            }else if(annotations$strand[i] == '+'){
                newRow[1,1] <- annotations$start[i]
                newRow[1,2] <- annotations$end[i]
                featureTable <- rbind(featureTable, newRow)
                if(identical(gene, getGeneName(annotations[i+1,])) 
                    && identical(annotations$strand[i], annotations$strand[i+1])){next}
            }
            qualifiers <- strsplit(annotations$attributes[i], ";")
            numQual <- length(qualifiers[[1]])
            qualifiers <- strsplit(qualifiers[[1]], "=")
            for(i in 1:numQual){
                if(any(str_to_lower(qualifiers[[i]][1]) == str_to_lower(illegalQualifiers))){next}
                addRow <- matrix('', 1, 5)
                addRow[1,4] <- qualifiers[[i]][1]
                addRow[1,5] <- qualifiers[[i]][2]
                featureTable <- rbind(featureTable, addRow)
            }
           
        }
        write.table(featureTable, paste(output, ".txt", sep=''), quote=FALSE, sep='\t', na='', row.names=FALSE, col.names=FALSE)
    }

    GB2FeatureTable <- function(gb, outPath, name){ 
       #gb <-readLines(paste(outFolderPath, genomes[i], ".gb", sep=''))
       #print('in MakeFeatureTable')
        keyFeatures <- c("assembly_gap", "C_region", "CDS", "centromere", "D-loop", "D_segment", "exon",
                        "gap", "gene", "iDNA", "intron", "J_segment", "mat_peptide", 
                        "misc_binding", "misc_difference", "misc_feature", "misc_recomb", "misc_RNA",
                        "misc_structure", "mobile_element", "modified_base", "mRNA", "ncRNA", "N_region",
                        "old_sequence", "operon", "oriT", "polyA_site", "precursor_RNA", "prim_transcript",
                        "primer_bind", "propeptide", "protein_bind", "regulatory", "repeat_region",
                        "rep_origin", "rRNA", "S_region", "sig_peptide", "source", "stem_loop", "STS",
                        "telomere", "tmRNA", "transit_peptide", "tRNA", "unsure", "V_region", "V_segment", 
                        "variation", "3'UTR", "5'UTR" )
        keyFeatures <- paste(" ", keyFeatures)
        features <- grep("FEATURES", gb) + 1
        endFeatures <- grep("ORIGIN", gb) -1  
        featureTable <- matrix('', (endFeatures - features), 5)
        j <- 2
        featureTable[1,1] <- paste(">Feature", name)

        for(i in features:endFeatures){
            if(j > nrow(featureTable)){
                addRow <- matrix('', 1, 5)
                featureTable <- rbind(featureTable, addRow)
            }
            if(stringr::str_detect(gb[i], "ORIGIN      				")){break}
    
            featureTable[j,1] <-  gb[i]
            line <- featureTable[j,1]
            keyFeatureCheck <- stringr::str_extract(line, "\\s\\s(\\w{3,15})")
            #Move key feature and start/ stop positions into proper place 
            if(any(str_to_lower(keyFeatureCheck) == str_to_lower(keyFeatures), na.rm = TRUE)){
                featureTable[j,3] <- stringr::str_extract(line, "(\\w{3,15})")
                if(grepl("misc_feature(\\s{4})", line)){ 
                    location <- stringr::str_extract(line, "([0-9]+)")
                    featureTable[j,2] <- location
                    featureTable[j,1] <- location
                }
                else{
                    while(stringr::str_detect(gb[i+1], "([0-9]+)(?=\\.\\..)")){
                        line <- paste(line, gb[i+1])
                        gb <- gb[-(i+1)]
                    }
                    first <- stringr::str_extract_all(line, "([0-9]+)(?=\\.\\..)")
                    second <- stringr::str_extract_all(line, "(?<=\\.\\.).([0-9]+)")
                    if(identical(first[[1]], character(0))){next}
                    for(k in 1:lengths(second)){
                        if(grepl("complement", line)){
                            featureTable[j, 2] <- first[[1]][k]
                            featureTable[j, 1] <- second[[1]][k]
                            if(k < lengths(second)) {j <- j+1}
                        }else{
                            featureTable[j,1] <- first[[1]][k]
                            featureTable[j,2] <- second[[1]][k]
                            if(k < lengths(second)) {j <- j+1}
                        }
                    }
                }
            }
            #Gets qualifier flags and breaks them up, putting into appropriate columns
            if(grepl("[[:blank:]]{21}\\/", featureTable[j,1])){
                featureTable[j,1] <- ''
                qualifier <- stringr::str_extract(line, "(?<=\\/)\\w{1,}(?==)")
                qualDescription <- stringr::str_extract(line, "(?<==).*") 
                if(is.na(qualifier)){next}
                featureTable[j,4] <- qualifier
                featureTable[j,5] <- qualDescription
            }
            #concatenates any multi line comments
            if(grepl("[[:blank:]]{21}\\w", featureTable[j,1]) == TRUE){
                line <- gsub("[[:blank:]]{21}", "", line)
                featureTable[j,1] <- ''
                featureTable[j-1,5] <- paste(featureTable[j-1,5], line, sep = '')
                next
            }
            j <- j +1 
        }
        #removes translation descriptors because those shouldn't be in feature tables
        for(i in 1:nrow(featureTable)){
            if(i > nrow(featureTable)){break}
            if(grepl("translation",featureTable[i,4]) |grepl("protein_id",featureTable[i,4])){
                featureTable <- featureTable[-i,]
            }
            if(all(is.na(featureTable[i,]))){
                featureTable <- featureTable[-i,]
            }
        }
        write.table(featureTable, file=paste(outPath, name, ".tbl", sep= ''), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
    }
    
    getNumberOfExons <- function(annotations, i, exons=0){
        #counts the number of exons in a particular gene
        for(j in (i[1]+1):(length(annotations$attributes)-1)){
            if(identical(annotations$attributes[j], annotations$attributes[i[1]]) && identical(annotations$strand[j],annotations$strand[i[1]])){
                i <- c(i, j)
            }
            else{
                return(i)
            }
        }
    }

    getCodingSequence <- function(annotations, dnachar, start, stop){
        #returns the sequence for a particular gene
        if(annotations[1,"strand"] == "-"){
            locus <- substr(dnachar, start, stop)
            locus <- DNAString(locus)
            codingSequence <- as.character(reverseComplement(locus))
        }
        if(annotations[1,"strand"] == "+"){
            codingSequence <- substr(dnachar, start, stop)
        }
        return(codingSequence)    
    }

    getGeneName <- function(annotation){
        attrib <- annotation[,'attributes']
        geneName <- str_extract(attrib, '(?<=Name=)(.[^;]*)') #stringr::str_extract(attrib, '(?<==).*(?=\\s)')
        return(geneName)
    }

    getGeneProduct <- function(annotation){
        attrib <- annotation[,'attributes']
        geneProduct <- str_extract(attrib, '(?<=product=)(.[^;]*)') #stringr::str_extract(attrib, '(?<==).*(?=\\s)')
        return(geneProduct)
    }

    codonGroup <- function(x){
        #takes a sequence and groups triplets into codons
        sequence.length <- nchar(x)
        triplet.starts <- seq(1, sequence.length, by=3) 
        lapply(triplet.starts, function(y){substring(x, y, y+2)})
        
    }

    CtoU_RNAediting <- function(codon){
        if(codon == 'ACG'){
            codon <- 'ATG'
        }
        if(codon == 'CAG'){
            codon <- 'TAG'
        }
        if(codon == 'CAA'){
            codon <- 'TAA'
        }
        if(codon == 'CGA'){
            codon <- 'TGA'
        }
        return(codon)
    }

    UtoC_RNAediting <- function(codon){
        if(codon == 'TAG'){
            codon <- 'CAG'
        }
        if(codon == 'TAA'){
            codon <- 'CAA'
        }
        if(codon == 'TGA'){
            codon <- 'CGA'
        }
        return(codon)
    }

    GFProteinFAsta <- function(output, genome, organism, gb){
        # >Seq1 [protein=neuropilin 1] [gene=Nrp1]
        # >ABCD [protein=merozoite surface protein 2] [gene=msp2] [protein_desc=MSP2]
        # >DNA.new [protein=breast and ovarian cancer susceptibility protein] [gene=BRCA1]
        protein <- as.data.frame(read.csv(paste(output, ".csv", sep='')))
        protein <- protein[-1,]
        protein <- arrange(protein, X)
        proteinFasta <- ""
        for(i in 1:nrow(protein)){
            line <- grep(protein$protein[i], gb)
            if(identical(line, integer(0))){
                cat("There seems to be a problem with the protein sequence of", as.character(protein$gene[i]), "please check its sequence in the protein fasta manually\n")
                next
            }
            endRange <- line[1]+10
            substrPositon <- grep('/product="', gb[(line[1]):endRange])
            gbLine <- line + substrPositon[1] - 1
            matches <- stringr::str_extract(gb[gbLine], '(?<=product=").*((?=\")|$)')
            if(!grepl('\"$', matches[1])){
                nextMatch <- stringr::str_extract(gb[gbLine+1], '(?<=\\s{21}).*(?=\")')
                matches <- paste(matches, nextMatch, sep = ' ')
            }else{
                matches<- stringr::str_extract(matches, '.*((?=\"))')
            }
            sequenceLength <- nchar(as.character(protein$protein[i]))
            gbRowLen <- seq(1, sequenceLength, by=70) 
            aaseq <- lapply(gbRowLen, function(y){substring(as.character(protein$protein[i]), y, y+69)})
            proteinFasta <-  paste(proteinFasta, ">lcl|", genome," [protein=", protein$proteinProduct[i], "] [gene=", protein$gene[i], "]","\n", paste(aaseq, collapse="\n"), "\n\n", sep = '')
        }
        writeLines(proteinFasta, paste(output, ".pep", sep=''))
    }

    GFFProteinFasta <- function(output, genome, organism){
        # >Seq1 [protein=neuropilin 1] [gene=Nrp1]
        # >ABCD [protein=merozoite surface protein 2] [gene=msp2] [protein_desc=MSP2]
        # >DNA.new [protein=breast and ovarian cancer susceptibility protein] [gene=BRCA1]
        protein <- as.data.frame(read.csv(paste(output, ".csv", sep='')))
        protein <- protein[-1,]
        protein <- arrange(protein, X)
        proteinFasta <- ""
        for(i in 1:nrow(protein)){
            sequenceLength <- nchar(as.character(protein$protein[i]))
            gbRowLen <- seq(1, sequenceLength, by=70) 
            aaseq <- lapply(gbRowLen, function(y){substring(as.character(protein$protein[i]), y, y+69)})
            proteinFasta <-  paste(proteinFasta, ">lcl|", genome," [protein=", protein$proteinProduct[i], "] [gene=", protein$gene[i], "]","\n", paste(aaseq, collapse="\n"), "\n\n", sep = '')
        }
        writeLines(proteinFasta, paste(output, ".pep", sep=''))
    }

    gb2fasta <- function(gb){
        #gets sequence data from GB file
        sequence <- ""
        for(i in 1:length(gb)){
            if(grepl("ORIGIN", gb[i])){
                for(j in i+1:length(gb)){
                    if(grepl("//", gb[j])){
                        return(sequence)
                    }else{
                        nucList <- unlist(stringr::str_extract_all(gb[j],'[a-z]'))
                        subSeq <- paste(nucList, collapse='')
                        sequence <- paste(sequence, subSeq, sep='')
                    }
                }
            }
        }
    }

    gff2fasta <- function(gffFile){
        sequence <- ""
        for(i in 1:length(gffFile)){
            if(grepl("#FASTA", gffFile[i])){
                startOfFas <- i   
                temp_gff <- gffFile[1:(i-1)]    
                for(j in i+2:length(gffFile)){
                    if(grepl("[A-z]", gffFile[j])){
                        nucList <- unlist(stringr::str_extract_all(gffFile[j],'[A-z]'))
                        subSeq <- paste(nucList, collapse='')
                        sequence <- paste(sequence, subSeq, sep='')
                    }else{
                        writeLines(temp_gff, "temp.gff")
                        return(sequence)
                    }
                }
            }
        }
    }

    seq2GFF <- function(output, dnachar, genome){
        sequence.length <- nchar(dnachar)
        seqByFourty <- seq(1, sequence.length, by=40) 
        splitUp <- lapply(seqByFourty[-length(seqByFourty)], function(y){substring(dnachar, y, y+39)})
        mashed <- paste(as.character(splitUp), sep="\n")
        write(paste("##FASTA\n>", genome, sep=''), file=paste(output, ".gff", sep=''),append=TRUE)
        write(mashed, file=paste(output, ".gff", sep=''),append=TRUE)
        write(substr(dnachar, (seqByFourty[length(seqByFourty)]), sequence.length),file=paste(output, ".gff", sep=''),append=TRUE )
    
    }

    restoreGFFheader <- function(gff, output){
        header <- gff[1]
        for(i in 2:length(gff)){
            if(grepl("##.*", gff[i])){
                header <- paste(header, gff[i], sep='\n')
            }else{
                newGFF <- readLines(paste(output, ".gff", sep=''))
                newGFF[1] <- paste(header, newGFF[1], sep='\n')
                writeLines(newGFF, paste(output, ".gff", sep=''))
                break
            }
        }
    }

    conventionalStops <- c('TAA','TAG','TGA')
    possibleEditedStops <- c('CAA','CAG','CGA')
    conventionalStarts <- c('ATG','GTG','TTG','ATT','CTG', 'ATC')
    possibleEditedStarts <- c('ACG') # 'ACC', 'CCG', 'GCG', 'CAG')

    for(i in 1:length(genomes)){
        #loops through the listed genome files and ReFerns each listed
        if(!is.null(gbFolderPath)){
            gff <- paste(gffFolderPath, genomes[i], ".gff", sep='')
            gb <- readLines(paste(gbFolderPath, genomes[i], ".gb", sep=''))
            gffFile <- readLines(gff)
            if(grepl("#FASTA", paste(gffFile, collapse=','))){
                sequence <- DNAString(gff2fasta(gffFile))
                dnachar <- as.character(sequence)
                annotations <- read.gff("temp.gff")
            }else{
                annotations <- read.gff(gff)
                if(is.null(fastaFolderPath)){
                    fasta <- readLines(paste(fastaFolderPath, genomes[i], ".fasta", sep=''))
                    sequence <- readDNAStringSet(fasta, format="fasta")
                }else{
                    sequence <- DNAString(gb2fasta(gb))
                }
            }
            dnachar <- as.character(sequence)
            output <- paste(outFolderPath, genomes[i], sep='')
            cat(genomes[i], "is being processed\n")
            editGB_File(genomes[i])
            gb <- readLines(paste(outFolderPath, genomes[i], ".gb", sep=''))
            GB2FeatureTable(gb, outFolderPath, genomes[i])
            organism <- stringr::str_extract(gb[1:20], "(?<=ORGANISM\\s\\s).*\\w")
            organism <- na.omit(organism)
            GFProteinFAsta(output, genomes[i], organism[1], gb)
            if(file.exists("temp.gff")){file.remove("temp.gff")}
        }else{
            gff <- paste(gffFolderPath, genomes[i], ".gff", sep='')
            if(is.null(fastaFolderPath)){
                gffFile <- readLines(gff)
                sequence <- DNAString(gff2fasta(gffFile))
                annotations <- read.gff("temp.gff")
            }else{
                fasta <- readLines(paste(fastaFolderPath, genomes[i], ".fasta", sep=''))
                sequence <- readDNAStringSet(fasta, format="fasta")
                annotations <- read.gff(gff)
            }
            gb <- NULL
            dnachar <- as.character(sequence)
            output <- paste(outFolderPath, genomes[i], sep='')
            cat(genomes[i], "is being processed\n")
            annotations <- edit_GFF_file(genomes[i])
            seq2GFF(output, dnachar, genomes[i])
            restoreGFFheader(gffFile, output)
            ##GET ORGANISM
            GFF2FeatureTable(annotations, output, genomes[i])
            #fix this
            GFFProteinFasta(output, genomes[i], genomes[i])
            if(file.exists("temp.gff")){file.remove("temp.gff")}
        }
        
    }
}
