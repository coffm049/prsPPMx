from ctypes import CDLL 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from xarray import DataArray
import pandas_plink as pdplink
import os
import argparse

try : 
  os.chdir("/home/christian/Research/Stat_gen/tools/MASH")
  from Simulate.simulation_helpers.Sim_generator import pheno_simulator
  from Simulate.summarizers.genotype_viz import plotClusters
  # from Simulate.summarizers.genotype_viz import plotClusters
  os.chdir("../../../Integrative/prsPPMx/")
except NameError :
  print("Already loaded appropriate modules")


def simulate_data(nSNPs=1000, nsubjects=500, nclusts=1, nphenos=2, shared=0.5, prop_causal=[0.25, 0.25], theta_alleles=[0.95, 0.25], h2Hom=0.8, h2Het=[0.1, 0.1],
                  prefix = "temp/simulation"):
   sim = pheno_simulator(nsubjects=nsubjects, nSNPs=nSNPs)
   sim.sim_sites(nsites = 1, nphenos= nphenos)
   sim.sim_pops(nclusts=nclusts, theta_alleles=theta_alleles, shared=shared)
   sim.sim_genos()
   sim.sim_pheno(h2Hom=h2Hom, h2Het=h2Het, prop_causal=prop_causal, alpha=-1, riskGroups = True, linearCombo = True)
   sim.save_plink(prefix = prefix)
   
rng = np.random.default_rng()


parser = argparse.ArgumentParser(description='Simulate data.')
parser.add_argument('--nSNPs', type=int, default=1000, help='Number of SNPs')
parser.add_argument('--nsubjects', type=int, default=500, help='Number of subjects')
parser.add_argument('--nclusts', type=int, default=1, help='Number of clusters')
parser.add_argument('--nphenos', type=int, default=2, help='Number of phenotypes')
parser.add_argument('--shared', type=float, default=0.5, help='Shared parameter')
parser.add_argument('--prop_causal', nargs='+', type=float, default=[0.25, 0.25], help='Proportion causal')
parser.add_argument('--theta_alleles', nargs='+', type=float, default=[0.95, 0.25], help='Theta alleles')
parser.add_argument('--h2Hom', type=float, default=0.8, help='Homogeneity')
parser.add_argument('--h2Het', nargs='+', type=float, default=[0.1, 0.1], help='Heterogeneity')
parser.add_argument('--prefix', type=str, default="temp/simulation", help='Prefix')

args = parser.parse_args()

simulate_data(nSNPs=args.nSNPs, nsubjects=args.nsubjects, nclusts=args.nclusts,
              nphenos=args.nphenos, shared=args.shared, prop_causal=args.prop_causal, theta_alleles=args.theta_alleles,
              h2Hom=args.h2Hom, h2Het=args.h2Het, prefix = "temp/simulation")