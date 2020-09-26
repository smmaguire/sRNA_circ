library(tidyverse)

bc1<-read_delim("barcode01_trimmed_processed.tab",delim="\t",
                 col_names = c("name","score","nedit","nmatch","nmismatch",
                               "trimmed_len","trimmed_seq")) %>%
  mutate(treatment = "no_amp")

bc2<-read_delim("barcode02_trimmed_processed.tab",delim="\t",
                col_names = c("name","score","nedit","nmatch","nmismatch",
                              "trimmed_len","trimmed_seq")) %>%
  mutate(treatment = "MDA")

bc3<-read_delim("barcode03_trimmed_processed.tab",delim="\t",
                col_names = c("name","score","nedit","nmatch","nmismatch",
                              "trimmed_len","trimmed_seq")) %>%
  mutate(treatment = "mouse")

nrow(filter(bc3,trimmed_len <= 23 & trimmed_len >= 18))/nrow(bc3)
full_data<-bind_rows(bc1,bc2)

ggplot(full_data,aes(x=treatment,y=nedit))+geom_violin()

nrow(filter(bc1, trimmed_len == 22))/nrow(bc1)
nrow(filter(bc2, nedit < 1))/nrow(bc2)

filter(bc1,nedit==0,trimmed_len != 22)

ggplot(full_data,aes(x=treatment,y=nmatch))+geom_violin()


ggplot(filter(full_data,nedit < 20),aes(x=treatment,y=nedit))+geom_violin()
ggplot(filter(full_data),aes(x=treatment,y=trimmed_len))+geom_violin() +
  ylab("Trimmed Length") + xlab(NULL)

ggplot(bc3,aes(x=trimmed_len))+geom_density()

ggplot(bc2,aes(x=trimmed_len))+geom_density()+geom_vline(xintercept=22)
median(bc3$trimmed_len)
