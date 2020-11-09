library(tidyverse)
library(readr)
library(cowplot)
data_loc<-"/Users/smaguire/Desktop/lab notebook - one note/small RNA nanopore/spade2/mirxplore_output"
files<-list.files(data_loc,full.names = T)

read_file<-function(file){
  barcode<-str_split(file,"barcode",simplify=T)[,2] %>%
    str_split("_covstats.txt",simplify=T)
  barcode<-as.numeric(barcode[,1])
  
  read_delim(file,col_names=c("ID","Avg_fold","Length","Ref_GC",
                              "Covered_percent","Covered_bases","Plus_reads",
                              "Minus_reads","Read_GC","Median_fold","std_dev"),
             delim="\t",
             col_types = "cdiddiiiddd") %>% 
    filter(ID != "#ID") %>%
    group_by(ID) %>%
    summarise(read_count = sum(Plus_reads,Minus_reads)) %>%
    ungroup() %>%
    mutate(total_reads=sum(read_count),
           expected_reads = total_reads/962,
           normalized_reads = read_count/expected_reads,
           barcode = barcode,
           amp = case_when(barcode == 4 | barcode == 5 ~ "No Amp",
                           barcode == 6 | barcode == 7 ~ "MDA"))
}

data<-map(files,read_file)

data<- bind_rows(data)

options(scipen = 999)
data %>%
  group_by(amp,ID) %>%
  summarise(mean_count = mean(normalized_reads)) %>%
  ungroup() %>%
  mutate(amp = factor(amp,levels=c("No Amp","MDA"))) %>%
  ggplot(aes(x=amp,y=mean_count,col=amp))+geom_jitter(size=0.4)+
  geom_boxplot(outlier.shape = NA,alpha=0.5,col="black")+
  geom_hline(yintercept=2,lty=2) +
  geom_hline(yintercept=0.5,lty=2) +
  geom_hline(yintercept=1,lty=2) + 
  scale_y_log10(breaks = c(0.001,0.01,.1,0.5,1,2,10,100)) +
  theme_cowplot() + 
  ylab("Mean Normalized Reads") +
  xlab(NULL) +
  guides(col="none")
  
p3<-data %>%
  group_by(amp,ID) %>%
  summarise(mean_count = mean(normalized_reads)) %>%
  ungroup() %>%
  pivot_wider(id_cols=ID,names_from=amp,values_from = mean_count) %>%
  ggplot(aes(x=`No Amp`,y=MDA))+geom_point() +
  scale_x_log10(limits=c(0.001,100))+
  scale_y_log10(limits=c(0.001,100)) +
  geom_smooth(method="lm") + 
  theme_cowplot()

p1<-data %>%
  filter(amp == "No Amp") %>%
  select(barcode,ID,normalized_reads) %>%
  pivot_wider(id_cols=ID,names_from=barcode,values_from = normalized_reads) %>%
  ggplot(aes(x=`4`,y=`5`))+geom_point() +
  scale_x_log10(limits=c(0.001,100))+
  scale_y_log10(limits=c(0.001,100)) + 
  geom_smooth(method="lm") +
  xlab("No Amp - Replicate 1") +
  ylab("No Amp - Replicate 2") + 
  theme_cowplot()

p2<-data %>%
  filter(amp == "MDA") %>%
  select(barcode,ID,normalized_reads) %>%
  pivot_wider(id_cols=ID,names_from=barcode,values_from = normalized_reads) %>%
  ggplot(aes(x=`6`,y=`7`))+geom_point() +
  scale_x_log10(limits=c(0.001,100))+
  scale_y_log10(limits=c(0.001,100))+
  geom_smooth(method="lm") +
  xlab("MDA - Replicate 1")+
  ylab("MDA - Replicate 2") + 
  theme_cowplot()

plot_grid(p1,p2,p3,align = "h")

mda.data<-
data %>%
  group_by(amp,ID) %>%
  summarise(mean_count = mean(normalized_reads)) %>%
  filter(amp == "MDA") %>%
  as.data.frame()


mda.data[(order(mda.data$mean_count,decreasing=T)),]


data %>%
  mutate(amp = factor(amp,levels=c("No Amp","MDA")),
         within2 = case_when(normalized_reads>=0.5 & normalized_reads <= 2 ~ TRUE,
                             TRUE ~ FALSE)) %>%
  group_by(amp,barcode,within2) %>%
  summarise(count=n()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(amp,barcode),names_from=within2,values_from=count) %>%
  mutate(total=`FALSE`+`TRUE`,
         percent_within2 = (`TRUE`/total)*100) %>%
  group_by(amp) %>%
  summarise(mean_percent_within2 = mean(percent_within2),
            sd_percent_within2 = sd(percent_within2)) %>%
  ggplot(aes(x=amp,y=mean_percent_within2,fill=amp))+
  geom_col()+
  geom_errorbar(aes(ymin=mean_percent_within2 - sd_percent_within2,
                    ymax=mean_percent_within2 + sd_percent_within2)) +
  theme_cowplot()+
  ylab("Percentage of miRNA within 2-fold of the expected value") +
  xlab(NULL) + ylim(c(0,80)) +
  guides(fill=F)
            