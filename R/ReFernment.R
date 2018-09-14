#' @title Annotate RNA editing in the GB files of plastomes
#'
#' @description annotates the GB files of plastomes with high levels of RNA
#' editing. Provides annotations to satisfy GenBank submission requirements. 
#' Takes as input the paths to the folders containing required file inputs
#' (a gff3 annotation and .gb file), the desired output path and the names of 
#' the files that the user wishes to annotate. Each filetype needn't be in 
#' seperate folders, but it is reccomended. This is especially true for output 
#' file path, because if it isn't sepereate from input gb files, the original 
#' files will be overwritten! Be sure that the file names for a respective 
#' genome are the same
#' across file types! 
#' 
#' @param gbFolderPath this is the path to the folder containing the gb file(s)
#' that you would like annotated. 
#' @param gffFolderPath this is the path to the folder containing the gff 
#' file(s) corresponding to the gff files you would like annotated. 
#' @param outFolderPath this is the path to the folder where the new genome
#' annotatons will be produced. The output is in gb format and the new 
#' annotations will be in the feature table. It is importatnt to note that you
#' should have the output file seperate from your input files, as this prevents
#' the script from overwriting the original files!
#' @param genomes this should be a list containing the filenames of all the 
#' genomes that you would like to annotate. Note that you should not include
#' the file endings of these files. 
#'
#' @return an annotated gb file. This file will have the required misc_feature
#' annotations as well as a translation corrected for RNA editing.   
#'
#' @author Tanner Robison
#' @examples
#' getwd()
#' genomes <- c("Asplenium_pek", "Woodwardia_uni")
#' gbFolder <- "/home/travis/build/TARobison/ReFernment/examples/GB/"
#' gffFolder <- "/home/travis/build/TARobison/ReFernment/examples/GFF/"
#' outputFolder <- "/home/travis/build/TARobison/ReFernment/examples/"
#' ReFernment(gbFolder, gffFolder, outputFolder, genomes)
#' 
#' @importFrom ape read.gff
#' @importFrom Biostrings DNAString DNAStringSet readDNAStringSet reverseComplement translate reverseComplement getGeneticCode
#' @importFrom stringr str_extract str_extract_all
#' @importFrom utils write.table
#' 
#' @export


ReFernment <- function(gbFolderPath, gffFolderPath, outFolderPath, genomes){
    
    editGB_File <- function(genomeName){
        #this function acts as a wrapper for many of the other functions contained within ReFernment. It adds RNA editing annotations as well as conceptual translations for the gb file. 
        #the protein fasta file is also made here. 
        proteinFasta <- ''
        for(i in 1:nrow(annotations)){
            gene <- getGeneName(annotations[i,])
            if(annotations$type[i] == "CDS"){
                exonIndex <- getNumberOfExons(annotations, i)
                if(identical(annotations$attributes[i], annotations$attributes[i-1]) && identical(annotations$strand[i], annotations$strand[i-1])){next}
                else if(identical(annotations$attributes[i], annotations$attributes[i+1]) && identical(annotations$strand[i], annotations$strand[i+1])){
                exons <- checkExons(annotations,i, dnachar,gb)  
                gb <- exons[[1]]
                protein <- exons[[2]]
                protein <- paste(">", genomeName, "|", gene, "\n", protein[1], "\n", sep = '')
                proteinFasta <- paste(proteinFasta, protein, sep='')
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

                    gb <- checkForInternalStops(CDS, annotations, i, gb)
                    
                    translation <- correctTranslation(annotations, i, dnachar, gb,CDS)
                    protein <- translation[1]
                    RNAediting <- translation[2]
                    gb <- annotateTranslation(protein, RNAediting, annotations,  gb,i)
                    protein <- paste(">", genomeName, "|", gene, "\n", protein, "\n", sep = '')
                    proteinFasta <- paste(proteinFasta, protein, sep='')
                    }
            }
        }
        writeLines(gb, paste(output, ".gb", sep=''))
        writeLines(proteinFasta, paste(output, "_protein.fasta", sep=''))
    }

    checkFrame <- function(annotations, i, CDS,gb){
        #I have no fucking clue what this does anymore
        gene <- getGeneName(annotations[i[1],])
        numExons <- length(i)
        if(numExons== 1){
            CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
        }   
        if(annotations$strand[i[numExons]] == '-'){
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
        }
        if(annotations$strand[i[numExons]] == '+'){
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
    }

    concatenateExons <- function(numExons, exonIndex, annotations, dnachar, i){
        #concatinates exons so we have a single coding sequence.
        firstExon <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
        for(j in 2:numExons){
                nextExon <- getCodingSequence(annotations[exonIndex[j],], dnachar, annotations$start[exonIndex[j]], annotations$end[exonIndex[j]])
                firstExon <- paste(firstExon, nextExon, sep = '')
        }
        return(firstExon)
    }

    checkExons <- function(annotations,i, dnachar, gb){
        #a wrapper function for genes with multiple invervals.
        exonIndex <- getNumberOfExons(annotations, i)
        numExons <- length(exonIndex)
        CDS <- concatenateExons(numExons, exonIndex, annotations, dnachar, i)

        if(nchar(CDS) %% 3 != 0){
            outCDS <- checkFrame(annotations, i, CDS,gb)
            gb <- outCDS[[1]]
            annotations <- outCDS[[2]]
            CDS <- concatenateExons(numExons, exonIndex, annotations, dnachar, i)
        }

        outStart <- checkStartCodon(CDS, annotations, exonIndex, gb)
        gb <- outStart[[1]]
        annotations <- outStart[[2]]
        CDS <- concatenateExons(numExons, exonIndex, annotations, dnachar, i)

        finalExon <- getCodingSequence(annotations[exonIndex[numExons],], dnachar, annotations$start[exonIndex[numExons]], annotations$end[exonIndex[numExons]])
        outStop <- checkStopCodon(finalExon, annotations,exonIndex[numExons], gb)
        gb <- outStop[[1]]
        annotations <- outStop[[2]]
        CDS <- concatenateExons(numExons, exonIndex, annotations, dnachar, i)

        gb <- checkForInternalStops(CDS, annotations,exonIndex,gb)
        translation <- correctTranslation(annotations, i, dnachar, gb,CDS)
        protein <- as.character(translation[1])
        RNAediting <- translation[2]
        gb <- annotateTranslation(protein, RNAediting, annotations,  gb,i)
        return(list(gb, protein))
    }

    correctTranslation <- function(annotations, i, dnachar, gb,CDS){
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

        annotatedSeq <- paste(codonList, collapse='')
    
        annotatedTranslation <- as.character(translate(DNAStringSet(annotatedSeq), getGeneticCode("11"),if.fuzzy.codon="solve"))
        
        out <- list(annotatedTranslation, RNAediting)

        return(out)
    }

    annotateTranslation <- function(annotatedTranslation, RNAediting, annotations,  gb,i){
        #edits GB file to add RNA editing qualifiers as well as conceptual translations.
        #browser()
        exonIndex <- getNumberOfExons(annotations, i)
        if(isTRUE(as.logical(RNAediting))){
            gbInsert <- paste('/translation="',annotatedTranslation, '"\n                     /exception="RNA editing"',sep='')
        }else{
            gbInsert <- paste('/translation="',annotatedTranslation, '"',sep='')
        }
        if(length(exonIndex) > 1){
            if(annotations[i[1], "strand"] == '-'){
            gbstring<- paste('(complement\\(',annotations[exonIndex[length(exonIndex)],"start"],'.*)',sep='')
            }
            if(annotations[i[1],"strand"] == '+'){
            gbstring<-paste('(CDS\\s*join\\(',annotations[i,"start"],'.*)',sep='')
            }
        }else{
            if(annotations[i[1], "strand"] == '-'){
                gbstring<- paste('(CDS\\s*complement\\(',annotations[i,"start"],'\\.\\.',annotations[i,"end"],'\\))',sep='')
            }
            if(annotations[i[1],"strand"] == '+'){
                gbstring<-paste('(CDS\\s*',annotations[i,"start"],'\\.\\.',annotations[i,"end"],')',sep='')
            }
        }
        newString <-  paste("\\1\n\                     ",gbInsert, sep='')
        gb <- sub(gbstring[1], newString[1], gb)
        
        return(gb)
    }

    checkStartCodon <- function(CDS, annotations,i,gb){
        #checks the first residue to see if it is a valid start codon. If not, checks to see whether RNA editing, or changing the gene boundaries restores start codon
        originalGB <- gb
        if(any(substr(CDS, 1, 3) == conventionalStarts)){
            out <- list(gb, annotations)
            return(out)
        }
        possibleEditedStarts <- c('ACG')
        if(any(substr(CDS, 1, 3) == possibleEditedStarts)){
            if(annotations[i[1], "strand"] == '+'){
                gb <- addFeature(annotations, i, annotations[i, "start"]+1, "CtoU",gb,"gene")
            }
            if(annotations[i[1], "strand"] == '-'){
                gb <- addFeature(annotations, i,annotations[i, "end"]-1, "CtoU",gb,"gene")
            }
            out <- list(gb, annotations)
            return(out)
        }else{
            gb <- shortenUpstream(annotations, dnachar, gb, CDS,i)
            if(identical(originalGB, gb[[1]])){
                gb <- extendUpstream(annotations, dnachar, gb[[1]], i)
            }
            return(gb)
        }   
    }

    checkStopCodon <- function(CDS, annotations,i, gb){
        #checks the final residue to see if it is a valid stop codon. If not, checks to see whether RNA editing, or changing the gene boundaries restores stop codon.
        exonIndex <- getNumberOfExons(annotations, i)
        numExons <- length(exonIndex)
        codonList <- codonGroup(CDS)
        originalGB <- gb
        if(any(substr(CDS, nchar(CDS)-2,nchar(CDS)) == (conventionalStops))){
            out <- list(gb, annotations)
            return(out) 
        } else if(any(substr(CDS, nchar(CDS)-2,nchar(CDS)) == possibleEditedStops)){
            if(annotations$strand[i[numExons]] == '+'){
                gb <- addFeature(annotations, i, annotations[i, "end"]-2, "CtoU",gb,"gene")
            }
            if(annotations$strand[i[numExons]] == '-'){
                gb <- addFeature(annotations, i,annotations[i, "start"]+2, "CtoU", gb, "gene")
            }
            out <- list(gb, annotations)
            return(out) 
        } else {
            gb <- shortenDownstream(annotations, dnachar, gb, CDS, i)
            if(identical(originalGB, gb[[1]])){
                gb <- extendDownstream(annotations, dnachar, gb[[1]], i)
            }
            return(gb)
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
                    gbstring <- paste(annotations$start[i[numExons]],'\\.\\.',sep='')
                    gb <- sub(gbstring, paste(newStop,'\\.\\.',sep=''), gb, perl = TRUE)
                    annotations$start[i[numExons]] <- newStop
                    out <- list(gb, annotations)
                    return(out)
                }
                if(annotations$strand[i[numExons]] == '+'){
                    dist <- j -(length(codonList) -10)
                    newStop <- annotations$end[i[numExons]] + dist*3
                    gbstring <- paste('\\.\\.',annotations$end[i[numExons]],sep='')
                    gb <- sub(gbstring, paste('\\.\\.', newStop,sep='') , gb, perl = TRUE)
                    annotations$end[i[numExons]] <- newStop
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
                    gbstring <- paste('\\.\\.',annotations$end[i[1]],sep='')
                    gb <- sub(gbstring, paste('\\.\\.', newStart[1],sep=''), gb, perl = TRUE)
                    annotations$end[i[1]] <- newStart
                    out <- list(gb, annotations)
                    return(out)
                }
                if(annotations$strand[i[1]] == '+'){
                    newStart <- annotations$start[i[1]] - count*3
                    gbstring <- paste(annotations$start[i[1]],'\\.\\.',sep='')
                    gb <- sub(gbstring, paste(newStart,'\\.\\.',sep='') , gb, perl = TRUE)
                    annotations$start[i[1]] <- newStart
                    out <- list(gb, annotations)
                    return(out)
                }  
            }
        }
        out <- list(gb, annotations)
        return(out)   
    }

    checkForInternalStops <- function(CDS, annotations ,i ,gb){
        #moves through each residue of a gene. If a stop codon is detected, function checks whether gene could be restored through RNA editing. 
        gene <- getGeneName(annotations[i,])
        codonList <-codonGroup(CDS)
        stops <- 0
        totalLength <- 0
        exonCount <- 1
            for(j in 1:(length(codonList)-1)){            
                if(annotations$strand[i[1]] == '-'){
                    if(annotations$end[i[exonCount]] - annotations$start[i[exonCount]] < j*3 && exonCount < length(i)){ 
                        totalLength <- annotations$end[i[exonCount]] - annotations$start[i[exonCount]] + totalLength +1
                        exonCount <- exonCount +1
                    }
                }else{
                    if(annotations$end[i[exonCount]] - annotations$start[i[exonCount]] < j*3 && exonCount < length(i)){
                        exonCount <- exonCount +1
                    }
                }
                if(any(codonList[j] == conventionalStops)){
                    if(annotations$strand[i[exonCount]] == '-'){
                        termStart <- annotations$end[i[exonCount]] - ((j-1)*3 - totalLength)
                        gb<-addFeature(annotations, i[1], termStart,"UtoC",gb,gene)
                        stops <- stops +1
                    }else{
                        termStart <- annotations$start[i[exonCount]] + ((j-1)*3 - totalLength)
                        gb<-addFeature(annotations, i[exonCount], termStart,"UtoC",gb,gene)
                        stops <- stops +1
                    }

                } 
            }
        if(stops > 5){
            gene <- getGeneName(annotations[i,])
            cat("There are a high number of edited Stops (", stops,") in", gene, "manually check to make sure frame is correct\n")
        }
        return(gb)
    }

    addFeature <- function(annotations,i,site,editType,gb,gene){
        #adds RNA editing qualifiers to gb file. 
        exonIndex <- getNumberOfExons(annotations,i)
        if(editType == 'CtoU'){
            note <- '/note="C to U RNA editing"'
        }
        if(editType == 'UtoC'){
            note <-  '/note="U to C RNA editing"'
        }
        newFeature <- paste('     misc_feature    ', site, '\n\                     ', note,'\n', sep ='')
        if(length(exonIndex) > 1){
            if(annotations[i[1], "strand"] == '-'){
                gbstring<- paste('(     CDS\\s*join\\(complement\\(',annotations[i,"start"],'.*)',sep='')
            }
            if(annotations[i[1],"strand"] == '+'){
                gbstring<-paste('(     CDS\\s*join\\(',annotations[i,"start"],'.*)',sep='')
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
        return(gb)
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
                    gbstring <- paste(annotations$start[i[numExons]],'\\.\\.',sep='')
                    gb <- sub(gbstring, paste(newStop,'\\.\\.',sep=''), gb, perl = TRUE)
                    annotations$start[i[numExons]] <- newStop
                    out <- list(gb, annotations)
                    return(out)
                }
                if(annotations$strand[i[numExons]] == '+'){
                    dist <- length(codonList) - j
                    newStop <- annotations$end[i[numExons]] - dist*3
                    gbstring <- paste('\\.\\.',annotations$end[i[numExons]],sep='')
                    gb <- sub(gbstring, paste('\\.\\.', newStop,sep='') , gb, perl = TRUE)
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
                    gbstring <- paste('\\.\\.',annotations[i,"end"],sep='')
                    gb <- sub(gbstring, paste('\\.\\.', newStart[1],sep=''), gb, perl = TRUE)
                    annotations[i,"end"] <- newStart
                    out <- list(gb, annotations)
                    return(out)
                }
                if(annotations[i[1],"strand"] == '+'){
                    newStart <- annotations[i,"start"] + (j-1)*3
                    gbstring <- paste(annotations[i,"start"], '\\.\\.',sep='')
                    gb <- sub(gbstring, paste(newStart, '\\.\\.',sep='') , gb, perl = TRUE)
                    annotations[i,"start"] <- newStart
                    out <- list(gb, annotations)
                    return(out)
                }    
            }
        }
        out <- list(gb, annotations)
        return(out)   
    }

    makeFeatureTable <- function(gb, outPath, name){ 
        keyFeatures <- c("assembly_gap", "C_region", "CDS", "centromere", "D-loop", "D_segment", 
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
            if(stringr::str_detect(gb[i], "ORIGIN      				")){break}
            featureTable[j,1] <-  gb[i]
            line <- featureTable[j,1]
            keyFeatureCheck <- stringr::str_extract(line, "\\s\\s(\\w{3,15})")
            #Move key feature and start/ stop positions into proper place 
            if(any(keyFeatureCheck == keyFeatures, na.rm = TRUE)){
                featureTable[j,3] <- stringr::str_extract(line, "(\\w{3,15})")
                if(grepl("misc_feature(\\s{4})", line)){ #
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
                qualDescription <- stringr::str_extract(line, "(?<==).*") #had a $ at end before
                if(is.na(qualifier)){next}
                featureTable[j,4] <- qualifier
                featureTable[j,5] <- qualDescription
            }
            #concatinates any multi line comments
            if(grepl("[[:blank:]]{21}\\w", featureTable[j,1]) == TRUE){
                line <- gsub("[[:blank:]]{21}", "", line)
                featureTable[j,1] <- ''
                featureTable[j-1,5] <- paste(featureTable[j-1,5], line, sep = '')
                next
            }
            j <- j +1 
        }
        #removes translation descriptors because those shouldn't be in feature tables, I guess
        for(i in 1:nrow(featureTable)){
            if(i > nrow(featureTable)){break}
            if(grepl("translation",featureTable[i,4]) |grepl("protein_id",featureTable[i,4])){
                featureTable <- featureTable[-i,]
            }
            if(all(is.na(featureTable[i,]))){
                featureTable <- featureTable[-i,]
            }
        }
        write.table(featureTable, file=paste(outPath, name, ".txt", sep= ''), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
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

    getCodingSequence <- function(annotations, dnachar, start,stop){
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
        geneName <- str_extract(attrib, '(?<==).*(?=\\s)') #stringr::str_extract(attrib, '(?<==).*(?=\\s)')
        return(geneName)
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

    conventionalStops <- c('TAA','TAG','TGA')
    possibleEditedStops <- c('CAA','CAG','CGA')
    conventionalStarts <- c('ATG','GTG','TTG','ATT','CTG')
    possibleEditedStarts <- c('ACG') # 'ACC', 'CCG', 'GCG', 'CAG')

    for(i in 1:length(genomes)){
        #loops through the listed genome files and ReFerns each listed
        gff <- paste(gffFolderPath, genomes[i], ".gff", sep='')
        gb <- readLines(paste(gbFolderPath, genomes[i], ".gb", sep=''))
        sequence <- DNAString(gb2fasta(gb))
        annotations <- read.gff(gff)
        dnachar <- as.character(sequence)
        output <- paste(outFolderPath, genomes[i], sep='')
        cat(genomes[i], "is being processed\n")
        editGB_File(genomes[i])
        gb <- readLines(paste(outFolderPath, genomes[i], ".gb", sep=''))
        makeFeatureTable(gb, outFolderPath, genomes[i])
    }
}

