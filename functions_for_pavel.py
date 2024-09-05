import numpy as np
import pandas as pd
import pickle
import subprocess
from Bio import SeqIO, AlignIO

### Building HMM model and aligning sequences to the model to build initial MSA for VAE training
def fix_seqs(input_string):
    seq_ = input_string.replace('.', '')
    seq_ = ''.join([char for char in seq_ if not char.islower()])
    return seq_

def hmmbuilding(input_msa, input_fasta, hmm_name, output_name,hmmer_path='/Users/dahala/Projects/HMMER/bin'):
    if not os.path.exists(f'/Users/dahala/GitHub/VAE-enzymes/hmm_model/{hmm_name}.hmm'):
        print('Building HMM model...')
        subprocess.run(f'{hmmer_path}/hmmbuild --amino /Users/dahala/GitHub/VAE-enzymes/hmm_model/{hmm_name}.hmm {input_msa} > /Users/dahala/GitHub/VAE-enzymes/hmm_model/{hmm_name}_hmmbuild.log', shell=True, executable="/bin/zsh")
    subprocess.run(f'cat {input_fasta} | {hmmer_path}/hmmalign --trim --outformat afa /Users/dahala/GitHub/VAE-enzymes/hmm_model/{hmm_name}.hmm - > /Users/dahala/GitHub/VAE-enzymes/vae_notebooks/datasets/{output_name}.afa', shell=True, executable="/bin/zsh")
    hmm_msa = AlignIO.read(f'/Users/dahala/GitHub/VAE-enzymes/vae_notebooks/datasets/{output_name}.afa', 'fasta')
    with open(f'/Users/dahala/GitHub/VAE-enzymes/vae_notebooks/datasets/{output_name}_fix.afa','w') as out_file:
        for prot in hmm_msa:
            out_file.write('>' + prot.id + '\n')
            out_file.write(fix_seqs(str(prot.seq)))
            out_file.write('\n')

hmmbuilding('/Users/dahala/GitHub/VAE-enzymes/datasets/GT1s_PtUGT1.fa', '/Users/dahala/GitHub/VAE-enzymes/datasets/GT1s_PtUGT1_raw.fas','EnzymeMiner_PtUGT1','EnzymeMiner_PtUGT1')

### Aligning sequences to the HMM model, filtering the MSA columns, and encoding the sequences to latent space

with open(run.pickles+ "/msa_columns.pkl", "rb") as input_file:
   msa_columns = pickle.load(input_file)
   
def fix_seqs(input_string):
    seq_ = input_string.replace('.', '')
    seq_ = ''.join([char for char in seq_ if not char.islower()])
    return seq_

def apply_msa_mask(msa_df, msa_columns):
    msa_df['trimmed_afa'] = msa_df['afa'].apply(lambda x: ''.join([x[i] for i in msa_columns]))
    return msa_df

def hmmer_align(fasta_name,hmm_model,out_path,hmmer_path='/Users/dahala/Projects/HMMER/bin',fasta_ext = 'faa'):
    subprocess.run(f'{hmmer_path}/hmmalign --trim --outformat afa /Users/dahala/GitHub/VAE-enzymes/hmm_model/{hmm_model}.hmm ../datasets/{fasta_name}.{fasta_ext} >> {out_path}/hmm/{fasta_name}.afa', shell=True, executable="/bin/zsh")
    return f'{out_path}/hmm/{fasta_name}.afa'

def aligned_df(afa_path,msa_columns=msa_columns):
    dict = {seq.description: fix_seqs(str(seq.seq)) for seq in SeqIO.parse(afa_path, "fasta")}
    df = pd.DataFrame.from_dict(dict, orient='index', columns=['afa'])
    df = apply_msa_mask(df, msa_columns)
    return df

def encode_custom_seqs(fasta_name,hmm_model,out_path,msa_columns=msa_columns,fasta_ext = 'faa',model= f'vae_fold_0.model',batch_size=1):
    afa_path = hmmer_align(fasta_name,hmm_model,out_path,fasta_ext = fasta_ext)
    df = aligned_df(afa_path,msa_columns)
    run.weights = f'{out_path}/model/{model}'
    latent_space = LatentSpace(run)
    mu1_ = []
    mu2_ = []
    for i in range(0, len(df), batch_size):
        latent_embeddings = latent_space.encode(df['trimmed_afa'].to_list()[i:i+batch_size])[0]
        mu1_.append(latent_embeddings[:,0])
        mu2_.append(latent_embeddings[:,1])
    df['mu1'] = np.concatenate(mu1_)
    df['mu2'] = np.concatenate(mu2_)
    return df

known_df = encode_custom_seqs('KnownGT1','EnzymeMiner_PtUGT1',run.results[:-7],fasta_ext='faa',model= f'vae_fold_0_e3999.model',batch_size=1)