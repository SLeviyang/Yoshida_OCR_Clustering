import pandas as pd
import sys

inputfile_name = sys.argv[1]

outputfile1 = sys.argv[2]
outputfile2 = sys.argv[3]

d = pd.read_csv(inputfile_name, sep='\t', header=None)
# shuffle
ds = d.sample(frac=1)

# split in two
nlines = ds.shape[0]
nsplit = round(nlines/2)
ds1 = ds.iloc[0:nsplit]
ds2 = ds.iloc[nsplit:nlines]

ds1.to_csv(outputfile1, index=False, header=None, sep="\t")
ds2.to_csv(outputfile2, index=False, header=None, sep="\t")
