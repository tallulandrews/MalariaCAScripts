# python3
import sys
import numpy as np
import pandas as pd 
import matplotlib
import matplotlib.pyplot as plt
import loompy
import velocyto as vcy
import logging

# save raw output 
vlm = vcy.VelocytoLoom("/lustre/scratch117/cellgen/team218/TA/Malaria_RNAVelocity/10XOutput/velocyto_output/possorted_genome_bam_UVJ4K.loom")

np.savetxt("out_S_matrix.csv", vlm.S, delimiter=",")
np.savetxt("out_U_matrix.csv", vlm.U, delimiter=",")
np.savetxt("out_A_matrix.csv", vlm.A, delimiter=",")

df = pd.DataFrame(vlm.ca)
df.to_csv("out_ca_matrix.csv")

df = pd.DataFrame(vlm.ra)
df.to_csv("out_ra_matrix.csv")

