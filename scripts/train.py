
import sys
import os
import os.path
import numpy as np
import pandas as pd
import torch
import pickle
import subprocess
from Bio import SeqIO
from importlib import reload

sys.path.append("/Users/dahala/GitHub/VAE-enzymes/vae_notebooks/")
sys.path.append("/Users/dahala/GitHub/VAE-enzymes/vae_notebooks/notebooks/")

#os.chdir("/Users/dahala/GitHub/VAE-enzymes/vae_notebooks/notebooks/")

import minimal_version.parser_handler
#os.chdir("/Users/dahala/GitHub/VAE-enzymes/vae_notebooks/")

from minimal_version.preprocess_msa import Preprocessor, weight_sequences
from minimal_version.msa import MSA
from minimal_version.utils import Capturing, store_to_pkl, store_to_fasta, load_from_pkl

from minimal_version.train import setup_train, Train
os.chdir("/Users/dahala/GitHub/VAE-enzymes/vae_notebooks/notebooks/")


#CONFIGURATION_FILE = "pfamGT1Small_filtered_hmm.json"
#CONFIGURATION_FILE = "msaEnzymeMiner_PtUGT1.json"
CONFIGURATION_FILE = "msaGASP_bigMSA.json"
run = minimal_version.parser_handler.RunSetup(CONFIGURATION_FILE)
print(f" Working with {CONFIGURATION_FILE} configuration file!")
PFAM_INPUT = False
'''
""" Configure model """
available_models = ["vae", "vae_conditional"]
model_selected = available_models[0]

""" Prepare data for training """
setup_train(run, model_selected)

"""
Train VAE with a preprocessed MSA
:param conf_path: path to the configuration file of the project
:return:
"""
# logging
train_log = open(os.path.join(run.logs, 'train_log.txt'), "w")
model = Train(run)
with Capturing() as output_train:
    model.train()
train_log.write(" Training \n" + "\n".join(output_train) + f'\n{"=" * 80}\n')
train_log.close()
'''