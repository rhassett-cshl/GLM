import pyreadr
import pandas as pd

# Specify the path to your RData file
rdata_file = "/home/rhassett/GLM/data/k562_epft_chr22_norm.RData"

# Read the RData file
result = pyreadr.read_r(rdata_file)

# The result is a dictionary-like object containing data frames and other objects
# You can access the data frames by their names
df = result["gb"]

print(df.head(5))