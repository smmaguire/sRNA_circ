library(Biostrings)
library(dplyr)
library(stringr)
# first argument will be the directory that contains the genbank files
# second argument will be the output directory
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
files<-list.files(pattern="repeat.gbk")

output_tb<-paste0(args[2],"/","temp_tb_",args[3],".tab")
output_fasta<-paste0(args[2],"/","temp_fa_",args[3],".fasta")

# function to rotate and align sequences to the adapter sequence
rotate_seq_align<-function(idf){
  seqs<-as.character(idf$unmasked_seq)
  file<-idf$file[1]
  seq.list<-c()
  names.list<-c()
  num<-1
  for(i in 1:length(seqs)){
    orig.seq<-seqs[i]
    seq<-seqs[i]
    seq_series<-letters[i]
    while(TRUE){
      seq.list<-c(seq.list,seq)
      names.list<-c(names.list,paste0(file,"_",seq_series,"_",num))
      new_seq<-str_sub(seq,start = 2)
      append<-str_sub(seq,end = 1)
      seq<-str_c(new_seq,append)
      num <- num+1
      if(seq == orig.seq){
        break
      }
    }
  }
  dna.seq<-DNAStringSet(seq.list)
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -2, baseOnly = TRUE)
  adapter<-DNAString("TACTGACCAGGACGACGACA")
  rc.adapter<-reverseComplement(adapter)
  alignment.a<-pairwiseAlignment(dna.seq, adapter, type = "global", substitutionMatrix = mat,
                                 gapOpening = 5, gapExtension = 0)
  alignment.b<-pairwiseAlignment(dna.seq, rc.adapter, type = "global", substitutionMatrix = mat,
                                 gapOpening = 5, gapExtension = 0)
  score.a<-score(alignment.a)
  score.b<-score(alignment.b)
 
  max.a<-max(score.a)
  max.b<-max(score.b)
 
  if(max.a > max.b){
    final.rotation<-seq.list[which(score.a == max.a)[1]]
    final.score<-max.a
  } else{
    final.rotation<-seq.list[which(score.b == max.b)[1]] %>%
    DNAString() %>%
    reverseComplement() %>%
    as.character()
    final.score<-max.b
  }
  idf<-mutate(idf,rotated_seq = final.rotation,
              aln_score = final.score)
  return(idf)
}

# main loop that goes through each genbank file, gets the info, rotates + aligns, outputs
# fasta and tab deliminated file.

for(gb in files){
  con2  <- file(gb, open = "r")
  full_length_line <- readLines(con2,n=1)
  split1<-str_split(str_squish(full_length_line)," ",simplify=T)
  full_seq_len<-as.numeric(split1[1,which(split1 == "bp")-1])
  read_id_line <- readLines(con2,n=1)
  split.read.id<-str_split(read_id_line,"DEFINITION",simplify = T)[2]
  read.id<-str_trim(split.read.id)
  close(con2)
  con  <- file(gb, open = "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(str_detect(line,"LOCUS")){
      split1<-str_split(str_squish(line)," ",simplify=T)
      seq_len<-as.numeric(split1[1,which(split1 == "bp")-1])
    }
    if(str_detect(line,"/period=")){
      split1<-str_split(str_split(line,"/period",simplify=T)[,2],"=",simplify = T)
      split2<-str_split(split1[,2],"\\\"",simplify=T)
      period<-as.numeric(split2[which(str_length(split2)>0)])
    }
    if(str_detect(line,"/rpt_num")){
      split1<-str_split(str_split(line,"/rpt_num",simplify=T)[,2],"=",simplify = T)
      split2<-str_split(split1[,2],"\\\"",simplify=T)
      rpt_num<-as.numeric(split2[which(str_length(split2)>0)])
    }
    if(str_detect(line,"/periodicity_score=")){
      split1<-str_split(str_split(line,"/periodicity_score",simplify=T)[,2],"=",simplify = T)
      split2<-str_split(split1[,2],"\\\"",simplify=T)
      period_score<-as.numeric(split2[which(str_length(split2)>0)])
    }
    if(str_detect(line,"/rpt_unit_seq=")){
      split1<-str_split(str_split(line,"/rpt_unit_seq",simplify=T)[,2],"=",simplify = T)
      split2<-str_split(split1[,2],"\\\"",simplify=T)
      rep_unit_seq<-split2[which(str_length(split2)>0)]
      while(TRUE){
        line<-readLines(con,n=1)
        if(str_detect(line,"/unmasked_rpt_unit_seq")){
          break
        }
        split3<-str_split(str_squish(line),"\\\"",simplify=T)[,1]
        rep_unit_seq<-str_c(rep_unit_seq,split3)
      }
    }
    if(str_detect(line,"/unmasked_rpt_unit_seq=")){
      split1<-str_split(str_split(line,"/unmasked_rpt_unit_seq",simplify=T)[,2],"=",simplify = T)
      split2<-str_split(split1[,2],"\\\"",simplify=T)
      unmasked_seq<-split2[which(str_length(split2)>0)]
      while(TRUE){
        line<-readLines(con,n=1)
        if(str_detect(line,"ORIGIN")){
          break
        }
        split3<-str_split(str_squish(line),"\\\"",simplify=T)[,1]
        unmasked_seq<-str_c(unmasked_seq,split3)
      }
    }
    if(str_detect(line,"ORIGIN")){
      break
    }
  }
  close(con)
  parsed.items<-data.frame(read.id,full_seq_len,seq_len,period,rpt_num,period_score,rep_unit_seq,unmasked_seq)
  parsed.items2<-rotate_seq_align(parsed.items)
  if(parsed.items2$aln_score <= 0){
    next
  }
  write.line<-paste(as.character(parsed.items2$read.id), as.character(parsed.items2$full_seq_len),
                    as.character(parsed.items2$seq_len), as.character(parsed.items2$period),
                    as.character(parsed.items2$rpt_num), as.character(parsed.items2$period_score),
                    as.character(parsed.items2$rep_unit_seq), as.character(parsed.items2$unmasked_seq),
                    as.character(parsed.items2$rotated_seq), as.character(parsed.items2$aln_score),
                    sep="\t")
  write.line<-paste(write.line,"\n")
  fasta_line<-c(">",as.character(parsed.items2$read.id),"\n",parsed.items2$rotated_seq,"\n")
  cat(write.line,file=output_tb,sep="\t",append=T)
  cat(fasta_line,file=output_fasta,sep="",append=T)
}
