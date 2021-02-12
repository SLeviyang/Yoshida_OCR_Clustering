library(varistran)
library(dplyr)

get_DATA_DIR <- function()
{
  return ("../../data_output/")
}


print("hello there")
infile <- paste(get_DATA_DIR(), "Yoshida/Yoshida_averaged_RNAseq.csv", sep="")
d = read.csv(infile, header=T, stringsAsFactors = F)
m = dplyr::select(d, -gene) %>% as.matrix()

y <- varistran::vst(m)
d_out <- cbind(d["gene"], as.matrix(y))
outfile <- paste(get_DATA_DIR(), "Yoshida/Yoshida_transformed_RNAseq.csv", sep="")
write.csv(d_out, outfile, row.names = F)