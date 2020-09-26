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
setwd(args[1])
fasta_file<-list.files(pattern="trimmed_output.fasta")

con <- file(fasta_file, open = "r")
barcode<-str_split(fasta_file,"_",simplify=T)[,2]
hsa_let7 <- DNAStringSet("TGAGGTAGTAGGTTGTATAGTT")
con2<-file(paste(barcode,"alignment_stats.txt",sep="_"),open = "w")

while ( TRUE ) {
  lines<-readLines(con, n = 2)
  if ( length(lines[1]) == 0 ) {
    break
  }
  name<-str_split(lines[1],">",simplify=T)[2]
  aln_global<-  pairwiseAlignment(hsa_let7, lines[2],
                                type = "global",
                                gapOpening = 1, gapExtension = 1)
  
  aln_local <- pairwiseAlignment(hsa_let7, lines[2],
                               type = "local",
                               gapOpening = 1, gapExtension = 1)
  
  aln_mouse <- pairwiseAlignment("AAGAAAGATTGCAAGAACTGCTAATTCATGCTTCCATGTTTAAAAACATGGCTTTCTTAC",
                                 lines[2],type="local",
                                 gapOpening = 1, gapExtension = 1)
  if(score(aln_local) < score(aln_mouse)){
    next
  }
  
  aln_global_stats<-aln_stats(aln_global)
  aln_local_stats<-aln_stats(aln_local)
  
  write_lines(paste(name,aln_global_stats,aln_local_stats,sep="\t"),path=con2,append=T)
}
close(con)
close(con2)
