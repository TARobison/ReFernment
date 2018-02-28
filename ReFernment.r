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
#checkExons
    #description needed
library(Biostrings)
library(genbankr)
library(annotate)
library(ape)




genomes <- c("AG-S2", "AG-S3","AG-S4","AG-S6","AG-S8","AG-S10","AG-S11","AG-S12","AG-S13","AG-S14","AG-S15","AG-S17", "AG-S18","AG-S20","AG-S21",
            "AG-S22","AG-S23","AG-S24","AG-S25","AG-S26","AG-S27","AG-S28","AG-S30","AG-S31","AG-S32","VittariaAppalachiana")

for(i in 1:length(genomes)){
    gff <- paste("C:\\Users\\tanner\\Documents\\RNAeditAnnotations\\GFF\\", genomes[i], ".gff", sep='') 
    fasta <- paste("C:\\Users\\tanner\\Documents\\RNAeditAnnotations\\FASTA\\", genomes[i], ".fasta", sep='')
    sequence <- readDNAStringSet(fasta, format="fasta")
    annotations <- read.gff(gff)
    dnachar <- as.character(sequence)
    gb <- readLines(paste("C:\\Users\\tanner\\Documents\\RNAeditAnnotations\\GB\\", genomes[i], ".gb", sep=''))
    output <- paste(      "C:\\Users\\tanner\\Documents\\RNAeditAnnotations\\Out\\", genomes[i], ".gb", sep='')
    print(genomes[i])
    editGB_File()
}

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
                gb <- checkExons(annotations,i, dnachar,gb)     
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

                gb <- checkForInternalStops(CDS, annotations,i,gb)
                }
        }
    }
    writeLines(gb, output)
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
    firstExon <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i])
    for(j in 2:numExons){
            nextExon <- getCodingSequence(annotations[exonIndex[j],], dnachar, annotations$start[exonIndex[j]], annotations$end[exonIndex[j]])
            firstExon <- paste(firstExon, nextExon, sep = '')
        }
    return(firstExon)
}


checkExons <- function(annotations,i, dnachar, gb){
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
    return(gb)
}

checkStartCodon <- function(CDS, annotations,i,gb){
    originalGB <- gb
    conventionalStarts <- c('ATG', 'GTG', 'TTG', 'ATT', 'CTG')
    if(any(substr(CDS, 1, 3) == conventionalStarts)){
        out <- list(gb, annotations)
        return(out)
    }
    possibleEditedStarts <- c('ACG', 'ACC', 'CCG', 'GCG', 'CAG')
    if(any(substr(CDS, 1, 3) == possibleEditedStarts)){
        if(annotations[i[1], "strand"] == '+'){
            gb <- addTranslationException(annotations, "MET", gb, annotations[i, "start"], annotations[i, "start"]+2,i)
        }
        if(annotations[i[1], "strand"] == '-'){
            gb <- addTranslationException(annotations, "MET", gb, annotations[i, "end"], annotations[i, "end"]-2,i)
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
    conventionalStops <- c('TAA','TAG','TGA')
    possibleEditedStops <- c('CAA', 'CAG', 'CGA')
    exonIndex <- getNumberOfExons(annotations, i)
    numExons <- length(exonIndex)
    codonList <- codonGroup(CDS)
    originalGB <- gb
    if(any(substr(CDS, nchar(CDS)-2,nchar(CDS)) == (conventionalStops))){
        out <- list(gb, annotations)
        return(out) 
    } else if(any(substr(CDS, nchar(CDS)-2,nchar(CDS)) == possibleEditedStops)){
        if(annotations$strand[i[numExons]] == '+'){
            gb <- addTranslationException(annotations, "TERM", gb, annotations[i, "end"], annotations[i, "end"]-2,i)
        }
        if(annotations$strand[i[numExons]] == '-'){
            gb <- addTranslationException(annotations, "TERM", gb, annotations[i, "start"]+2, annotations[i, "start"],i)
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
    if(annotations$strand[i[1]] == '-'){
        CDS <- getCodingSequence(annotations[i,], dnachar, (annotations$start[i]-30), annotations$end[i])
    }
    if(annotations$strand[i[1]] == '+'){
        CDS <- getCodingSequence(annotations, dnachar, annotations$start[i], annotations$end[i] +30)
    }
    codonList <- codonGroup(CDS)
    conventionalStops  <-  c('TAA','TAG','TGA')
    possibleEditedStops <- c('CAA','CAG','CGA')
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
    if(annotations$strand[i[1]] == '-'){
        CDS <- getCodingSequence(annotations[i,], dnachar, annotations$start[i], annotations$end[i]+18)
    }
    if(annotations$strand[i[1]] == '+'){
        CDS <- getCodingSequence(annotations, dnachar, annotations$start[i] -18, annotations$end[i])
    }
    codonList <- codonGroup(CDS)
    conventionalStarts  <-  c('ATG', 'GTG', 'TTG', 'ATT', 'CTG')
    possibleEditedStarts <- c('ACG', 'ACC', 'CCG', 'GCG', 'CAG')
    count <- 0
    for(j in 6:1){
        count <- count+1
         if(any(codonList[j] == conventionalStarts)){
             if(annotations$strand[i[1]] == '-'){
                newStart <- annotations$end[i[1]] + count*3
                gbstring <- paste('\\.\\.',annotations$end[i[1]],sep='')
                gb <- sub(gbstring, paste('\\.\\.', newStart,sep=''), gb, perl = TRUE)
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
        cat("There are a high number of edited Stops (", stops,") in", gene, "manually check to make sure frame is correct\n")
    }
    return(gb)
}

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
    codonList <-codonGroup(CDS)
    conventionalStarts <- c('ATG', 'GTG', 'TTG', 'ATT', 'CTG')
    for(j in 1:5){
        if(any(codonList[j] == conventionalStarts)){
            if(annotations[i[1],"strand"] == '-'){
                newStart <- annotations[i,"end"] - (j-1)*3
                gbstring <- paste('\\.\\.',annotations[i,"end"],sep='')
                gb <- sub(gbstring, paste('\\.\\.', newStart,sep=''), gb, perl = TRUE)
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

getNumberOfExons <- function(annotations, i, exons=0){
    for(j in i+1:length(annotations$attributes)){
        if(identical(annotations$attributes[j], annotations$attributes[i[1]]) && identical(annotations$strand[j],annotations$strand[i[1]])){
            i <- c(i, j)
        }
        else{
            return(i)
        }
    }
   # return(i)
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