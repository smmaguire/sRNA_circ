library(Biostrings)
library(dplyr)
library(stringr)
library(readr)

aln_stats<-function(alignment,reference){
  ref_len<-str_length(reference)
  edits<-nedit(alignment)
  matches<-nmatch(alignment)
  mismatches<-nmismatch(alignment)
  per_id<-round(((ref_len-edits)/ref_len)*100,2)
  aln_score<-round(score(alignment),2)
  paste(matches,mismatches,edits,aln_score,per_id,sep="\t")
}

args <- commandArgs(trailingOnly = TRUE)
file<-args[1]

start<-F
seq<-""
hsa_let7 <- DNAStringSet("TGAGGTAGTAGGTTGTATAGTT")

con <- file(file, open = "r")
write_file <- paste(args[2],"_alignment_stats.txt",sep="_")
pb <- txtProgressBar(min=0,max=as.numeric(args[3]),style=3)
line_count<-0
while ( TRUE ) {
  line = readLines(con, n = 1)
  line_count<-line_count+1
  setTxtProgressBar(pb,line_count)
  if ( length(line) == 0 ) {
     break
   }
  if(start & str_detect(line,">")){
   aln_global<-  pairwiseAlignment(hsa_let7, seq,
                                type = "global",
                                gapOpening = 1, gapExtension = 1)
  
  aln_local <- pairwiseAlignment(hsa_let7, seq,
                               type = "local",
                               gapOpening = 1, gapExtension = 1)
  
  aln_mouse <- pairwiseAlignment("AAGAAAGATTGCAAGAACTGCTAATTCATGCTTCCATGTTTAAAAACATGGCTTTCTTAC",
                                 seq,type="local",
                                 gapOpening = 1, gapExtension = 1)
  
  if(score(aln_local) < score(aln_mouse)){
    next
  }
  
  aln_global_stats<-aln_stats(aln_global,hsa_let7)
  aln_local_stats<-aln_stats(aln_local,hsa_let7)
  write.line<-paste(name,aln_global_stats,aln_local_stats,sep="\t")
  write.line<-paste0(write.line,"\n")
  cat(write.line,file=write_file,append = T)
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
close(pb)
close(con)
