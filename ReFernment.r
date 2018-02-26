#check start codon
    #This function checks the start codon, if it is a valid conventional start codon,
    #the function is exited. If the first residues become valid through RNA editing then
    #the annotations are altered accordingly. If it isn't either of the above cases, you will have to
    #manually edit the CDS
#check stop codon
    #description needed
#checkInternalStops
    #description needed
#getCodingSequence
    #this function uses the start and stop positions of a particular cds to extract the
    #cds, and if it is on the minus strand, it will get the compliment of the sequence, 
    #so that the correct translation can be given
#concatenateExons
    #description needed
library(Biostrings)
library(genbankr)
library(annotate)
library(ape)

# gff <- "C:\\Users\\tanner\\Documents\\RNAeditAnnotations\\GFF\\AG-S3.gff" 
# fasta <- "C:\\Users\\tanner\\Documents\\RNAeditAnnotations\\FASTA\\AG-S3.fasta"
# sequence <- readDNAStringSet(fasta, format="fasta")
# annotations <- read.gff(gff)
# dnachar <- as.character(sequence)
# gb <- readLines("C:\\Users\\tanner\\Documents\\RNAeditAnnotations\\GB\\AG-S3.gb")

# args = commandArgs(trailingOnly=TRUE)

# gff <- args[1]
# annotations <- read.gff(gff)
# fasta <- args[2]
# sequence <- readDNAStringSet(fasta, format="fasta")
# dnachar <- as.character(sequence)
# gb <- readLines(args[3])

editGB_File <- function(){

    for(i in 1:nrow(annotations)){
        gene <- getGeneName(annotations[i,])
        if(annotations$type[i] == "CDS"){
            exonIndex <- getNumberOfExons(annotations, i)
            if(identical(annotations$attributes[i], annotations$attributes[i-1]) && identical(annotations$strand[i], annotations$strand[i-1])){
                next        
            }else if(identical(annotations$attributes[i], annotations$attributes[i+1]) && identical(annotations$strand[i], annotations$strand[i+1])){
                gb <- concatenateExons(annotations,i, dnachar,gb)       
            }else{
                CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
                if(nchar(CDS) %% 3 != 0){
                    out <- checkFrame(annotations, i, CDS,gb)
                    gb <- out[[1]]
                    if(annotations$strand[i] == '-'){
                        annotations$start[i] <- out[[2]]
                    }else{
                        annotations$end[i] <- out[[2]]
                    }
                    CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
                }
                
                gb <- checkStartCodon(CDS, annotations,i,gb)
                if(length(gb) < 4){
                    if(annotations$strand[i] == '-'){
                        annotations$end[i] <- gb[[2]]
                    }else{
                        annotations$start[i] <- gb[[2]]
                    }
                    gb <- gb[[1]]
                    CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
                }

                out <- checkStopCodon(CDS, annotations,i, gb)
              
                if(length(out) < 4){
                    gb <- out[[1]]
                    if(annotations$strand[i] == '-'){
                        annotations$start[i] <- out[[2]]
                    }else{
                        annotations$end[i] <- out[[2]]
                    }
                    CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
                }else{
                    gb <- out
                }
                gb <- checkForInternalStops(CDS, annotations,i,gb)
                }
        }
    }
    writeLines(gb, "C:\\Users\\tanner\\Documents\\RNAeditAnnotations\\testyBoi.gb")#con = args[4])
}

checkFrame <- function(annotations, i, CDS,gb){
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
        out <- list(gb,annotations$start[i[numExons]])
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
        out <- list(gb, annotations$end[i[numExons]])
        return(out)
    }
}

concatenateExons <- function(annotations,i, dnachar, gb){
    exonIndex <- getNumberOfExons(annotations, i)
    numExons <- length(exonIndex)
    if(getGeneName(annotations[i,]) == "ndhA"){
    }
    if(annotations$strand[i[1]] == "-"){
        firstExon <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
        out <- checkStartCodon(firstExon, annotations, exonIndex, gb)
        if(length(out) < 4){
            annotations$start[exonIndex[numExons]] <- as.numeric(out[[2]])
            gb <- out[[1]]
        }else{
            gb <- out
        }
        for(j in 2:numExons){
            nextExon <- getCodingSequence(annotations[exonIndex[j],], dnachar, annotations$start[exonIndex[j]], annotations$end[exonIndex[j]])
            firstExon <- paste(firstExon, nextExon, sep = '')
        }

        out <- checkFrame(annotations, exonIndex, firstExon,gb)
        annotations$start[exonIndex[numExons]] <- as.numeric(out[[2]])
        gb <- out[[1]]
        
        #I think the problem is that this function is too damn big. It needs to be refactored into a few smaller functions. 
        out <- checkStopCodon(nextExon, annotations,exonIndex[numExons], gb)
        if(length(out) < 4){
            annotations$start[exonIndex[numExons]] <- as.numeric(out[[2]])
            gb <- out[[1]]
        }else{
            gb <- out
        }
        gb <- checkForInternalStops(firstExon, annotations,exonIndex,gb)
        return(gb)
    }else{
        firstExon <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
        out <- checkStartCodon(firstExon, annotations,i, gb)
        if(length(out) < 4){
            annotations$start[exonIndex[numExons]] <- as.numeric(out[[2]])
            gb <- out[[1]]
        }else{
            gb <- out
        }
        for(j in 2:numExons){
            nextExon <- getCodingSequence(annotations[exonIndex[j],], dnachar, annotations$start[exonIndex[j]], annotations$end[exonIndex[j]])
            firstExon <- paste(firstExon, nextExon, sep='')
        }

        out <- checkFrame(annotations, exonIndex, firstExon,gb)
        annotations$end[exonIndex[numExons]] <- as.numeric(out[[2]])
        gb <- out[[1]]

        out <- checkStopCodon(nextExon, annotations,j, gb)
        if(length(out) < 4){
            annotations$end[exonIndex[numExons]] <- as.numeric(out[[2]])
            gb <- out[[1]]
        }else{
            gb <- out
        }
        gb <- checkForInternalStops(firstExon, annotations,i, gb)
        return(gb)
        }
    return(gb)
}

checkStartCodon <- function(CDS, annotations,i,gb){
    originalGB <- gb
    conventionalStarts <- c('ATG', 'GTG', 'TTG', 'ATT', 'CTG')
    if(any(substr(CDS, 1, 3) == conventionalStarts)){
        return(gb)
    }
    possibleEditedStarts <- c('ACG', 'ACC', 'CCG', 'GCG', 'CAG')
    if(any(substr(CDS, 1, 3) == possibleEditedStarts)){
        if(annotations[i[1], "strand"] == '+'){
            gb <- addTranslationException(annotations, "MET", gb, annotations[i, "start"], annotations[i, "start"]+2,i)
        }
        if(annotations[i[1], "strand"] == '-'){
            gb <- addTranslationException(annotations, "MET", gb, annotations[i, "end"], annotations[i, "end"]-2,i)
        }
        return(gb)
    }else{
        gb <- shortenUpstream(annotations, dnachar, gb, CDS,i)
         if(identical(originalGB, gb)){
            gb <- extendUpstream(annotations, dnachar, gb, i)
        }
    } 
    return(gb)   
}

checkStopCodon <- function(CDS, annotations,i, gb){
    conventionalStops <- c('TAA','TAG','TGA')
    possibleEditedStops <- c('CAA', 'CAG', 'CGA')
    exonIndex <- getNumberOfExons(annotations, i)
    numExons <- length(exonIndex)
    codonList <- codonGroup(CDS)
    originalGB <- gb
    if(any(substr(CDS, nchar(CDS)-2,nchar(CDS)) == (conventionalStops))){
        return(gb)
    } else if(any(substr(CDS, nchar(CDS)-2,nchar(CDS)) == possibleEditedStops)){
        if(annotations$strand[i[numExons]] == '+'){
            gb <- addTranslationException(annotations, "TERM", gb, annotations[i, "end"], annotations[i, "end"]-2,i)
        }
        if(annotations$strand[i[numExons]] == '-'){
            gb <- addTranslationException(annotations, "TERM", gb, annotations[i, "start"]+2, annotations[i, "start"],i)
        }
        return(gb)
    } else {
        gb <- shortenDownstream(annotations, dnachar, gb, CDS, i)
        if(identical(originalGB, gb)){
            gb <- extendDownstream(annotations, dnachar, gb, i)
        }
        return(gb)
    }
}

extendDownstream <- function(annotations, dnachar, gb, i){
    if(annotations$strand[i[1]] == '-'){
        CDS <- getCodingSequence(annotations[i,], dnachar, (annotations$start[i]-15), annotations$end[i])
    }
    if(annotations$strand[i[1]] == '+'){
        CDS <- getCodingSequence(annotations, dnachar, annotations$start[i], annotations$end[i] +15)
    }
    codonList <- codonGroup(CDS)
    conventionalStops  <-  c('TAA','TAG','TGA')
    possibleEditedStops <- c('CAA','CAG','CGA')
    exonIndex <- getNumberOfExons(annotations, i)
    numExons <- length(exonIndex)
    for(j in (length(codonList)-5):((length(codonList)))){
        if(any(codonList[j] == conventionalStops)){
            if(annotations$strand[i[numExons]] == '-'){
                dist <- j - (length(codonList)-5)
                newStop <- annotations$start[i[numExons]] - dist*3
                gbstring <- paste(annotations$start[i[numExons]],'\\.\\.',sep='')
                gb <- sub(gbstring, paste(newStop,'\\.\\.',sep=''), gb, perl = TRUE)
                annotations$start[i[numExons]] <- newStop
                return(list(gb, newStop))
            }
            if(annotations$strand[i[numExons]] == '+'){
                   
                dist <- j -(length(codonList) -5)
                newStop <- annotations$end[i[numExons]] + dist*3
                gbstring <- paste('\\.\\.',annotations$end[i[numExons]],sep='')
                gb <- sub(gbstring, paste('\\.\\.', newStop,sep='') , gb, perl = TRUE)
                return(list(gb, newStop))
            }  
        }
    }
    return(gb)   
}

extendUpstream <- function(annotations,dnachar, gb, i){
    if(annotations$strand[i[1]] == '-'){
        CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i]+15)
    }
    if(annotations$strand[i[1]] == '+'){
        CDS <- getCodingSequence(annotations, dnachar, annotations$start[i] -15, annotations$end[i])
    }
    codonList <- codonGroup(CDS)
    conventionalStarts  <-  c('ATG', 'GTG', 'TTG', 'ATT', 'CTG')
    possibleEditedStarts <- c('ACG', 'ACC', 'CCG', 'GCG', 'CAG')
    count <- 0
    for(j in 5:1){
        count <- count+1
         if(any(codonList[j] == conventionalStarts)){
             if(annotations$strand[i[1]] == '-'){
                newStart <- annotations$end[i[1]] + count*3
                gbstring <- paste('\\.\\.',annotations$end[i[1]],sep='')
                gb <- sub(gbstring, paste('\\.\\.', newStart,sep=''), gb, perl = TRUE)
                return(list(gb, newStart))
            }
            if(annotations$strand[i[1]] == '+'){
                newStart <- annotations$start[i[1]] - count*3
                gbstring <- paste(annotations$start[i[1]],'\\.\\.',sep='')
                gb <- sub(gbstring, paste(newStart,'\\.\\.',sep='') , gb, perl = TRUE)
                return(list(gb, newStart))
            }  
        }
    }
    return(gb)   
}

checkForInternalStops <- function(CDS, annotations ,i ,gb){
    gene <- getGeneName(annotations[i,])
    codonList <-codonGroup(CDS)
    conventionalStops <- c('TAA','TAG','TGA')
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
                    termStop  <- annotations$end[i[exonCount]] - ((j-1)*3 - totalLength) - 2 
                    stops <- stops +1
                }else{
                    termStart <- annotations$start[i[exonCount]] + ((j-1)*3 - totalLength)
                    termStop  <- annotations$start[i[exonCount]] + ((j-1)*3 - totalLength) + 2 
                    stops <- stops +1
                }
                gb <- addTranslationException(annotations,"other", gb, termStart, termStop,i)
            } 
        }
    if(stops > 5){
        gene <- getGeneName(annotations[i,])
        cat("There are a high number of edited Stops in", gene, "manually check to make sure frame is correct\n")
    }
    return(gb)
}

#changed this from 5 codons to 7. If you notice anything wonky, it might be from this
shortenDownstream <- function(annotations, dnachar, gb, CDS,i){
    codonList <-codonGroup(CDS)
    conventionalStops <- c('TAA','TAG','TGA')
    exonIndex <- getNumberOfExons(annotations, i)
    numExons <- length(exonIndex)
    for(j in (length(codonList)-7):((length(codonList)-1))){
        if(any(codonList[j] == conventionalStops)){
            if(annotations$strand[i[numExons]] == '-'){
                dist <- length(codonList) - j
                newStop <- annotations$start[i[numExons]] + dist*3
                gbstring <- paste(annotations$start[i[numExons]],'\\.\\.',sep='')
                gb <- sub(gbstring, paste(newStop,'\\.\\.',sep=''), gb, perl = TRUE)
                ll <- list(gb,newStop)
                return(ll)
            }
            if(annotations$strand[i[numExons]] == '+'){
                dist <- length(codonList) - j
                newStop <- annotations$end[i[numExons]] - dist*3
                gbstring <- paste('\\.\\.',annotations$end[i[numExons]],sep='')
                gb <- sub(gbstring, paste('\\.\\.', newStop,sep='') , gb, perl = TRUE)
                return(list(gb,newStop))
            }  
        }
    }
    return(gb)   
}

shortenUpstream <- function(annotations, dnachar, gb, CDS,i){
    codonList <-codonGroup(CDS)
    conventionalStarts <- c('ATG', 'GTG', 'TTG', 'ATT', 'CTG')
    for(j in 1:5){
        if(any(codonList[j] == conventionalStarts)){
            if(annotations[i[1],"strand"] == '-'){
                newStart <- annotations[i,"end"] - (j-1)*3
                gbstring <- paste('\\.\\.',annotations[i,"end"],sep='')
                gb <- sub(gbstring, paste('\\.\\.', newStart,sep=''), gb, perl = TRUE)
                return(list(gb, newStart))
            }
            if(annotations[i[1],"strand"] == '+'){
                newStart <- annotations[i,"start"] + (j-1)*3
                gbstring <- paste(annotations[i,"start"], '\\.\\.',sep='')
                gb <- sub(gbstring, paste(newStart, '\\.\\.',sep='') , gb, perl = TRUE)
                return(list(gb, newStart))
            }    
        }
    }
    return(gb)   
}

getNumberOfExons <- function(annotations, i, exons=0){
    for(j in i+1:length(annotations$attributes)){
        if(identical(annotations$attributes[j], annotations$attributes[i[1]]) && identical(annotations$strand[j],annotations$strand[i[1]])){
            i <- c(i, j)
        }
    }
    return(i)
}

getCodingSequence <- function(annotations, dnachar, start,stop){
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
    geneName <- stringr::str_extract(attrib, '(?<==).*(?=\\s)')
    return(geneName)
}

codonGroup <- function(x){
    sequence.length <- nchar(x)
    triplet.starts <- seq(1, sequence.length, by=3) 
    lapply(triplet.starts, function(y){substring(x, y, y+2)})
    
}

addTranslationException <- function(annotations, aa, gb, start, stop,i){
    exceptSite <- paste('/transl_except="(pos:',start,'-',stop,',aa:',aa,')"',sep='' )
    if(length(i) > 1){
        i <- i[length(i)]
        gbstring<- paste('(',annotations$start[i],'\\.\\.',annotations$end[i],'\\)\\))',sep='')
        }else{
        if(annotations[i[1], "strand"] == '-'){
            gbstring<- paste('(CDS\\s*complement\\(',annotations[i,"start"],'\\.\\.',annotations[i,"end"],'\\))',sep='')
        }
        if(annotations[i[1],"strand"] == '+'){
            gbstring<-paste('(CDS\\s*',annotations[i,"start"],'\\.\\.',annotations[i,"end"],')',sep='')
        }
    }
    newString <-  paste("\\1\n\                     ",exceptSite,sep='')
    gene <- getGeneName(annotations[i,])
    gb <- sub(gbstring, newString[1], gb)
    return(gb)
}

editGB_File()

