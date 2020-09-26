library(Biostrings)
library(dplyr)
library(stringr)
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
files<-list.files(pattern="trimmed")

for(file in files){
  con <- file(file, open = "r")
  barcode<-str_split(file,"_",simplify=T)[,2]
  start<-F
  seq<-""
  hsa_let7 <- DNAStringSet("TGAGGTAGTAGGTTGTATAGTT")
  rc_hsa_let7 <- reverseComplement(hsa_let7)
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
  while ( TRUE ) {
    line = readLines(con, n = 1)
    line
    if ( length(line) == 0 ) {
      break
    }
    if(start & str_detect(line,">")){
      aln_fwd <-  pairwiseAlignment(hsa_let7, seq,
                                  type = "global", substitutionMatrix = mat,
                                  gapOpening = 1, gapExtension = 1)
 
    aln_rev <- pairwiseAlignment(rc_hsa_let7, seq,
                                 type = "global", substitutionMatrix = mat,
                                 gapOpening = 1, gapExtension = 1)
    if(score(aln_fwd) > score(aln_rev)) {
      aln_go <- aln_fwd
    } else{
      aln_go <- aln_rev
    }
    write.line<-paste(name,score(aln_go),nedit(aln_go), nmatch(aln_go), nmismatch(aln_go),
          nchar(seq), seq, sep = "\t")
    write.line<-paste0(write.line,"\n")
    cat(write.line,file=paste0(barcode,"_","processed.tab"),append = T)
    seq<-""
    start<-!start
    }
    if(!start & str_detect(line,">")){
      name<-str_split(line,">",simplify=T)[2]
      start=!start
    }
    if(start & !str_detect(line,">")){
      seq<-paste0(seq,line)
    }
}
close(con)
}
