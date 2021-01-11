library(dplyr)

# read in chromVar mouse v2 cisDB motifs used by Yoshida et al as TF motifs
# downloaded from
# https://github.com/buenrostrolab/chromVARmotifs/blob/master/data/mouse_pwms_v2.rda
load("../../data_input/mouse_pwms_v2.rda")
motifs <- mouse_pwms_v2
motif_names <- sapply(motifs, function(m) name(m)) %>% setNames(NULL)
names(motifs) <- motif_names

# Yoshida et al restrict to these 430 TF 
d <- read.csv("../../data_input/Yoshida_Table_S5.csv", header=T, stringsAsFactors = F)
Yoshida_motifs <- motifs[is.element(motif_names, d$TF)]

save(Yoshida_motifs, file="../../data_input/Yoshida_TF_PWM.rda")