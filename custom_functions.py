import numpy as np
import pandas as pd
import pickle
import subprocess
import os
from Bio import SeqIO, AlignIO

class HMMmodule:

    def __init__(self,run_config):

        self.muscle_path = run_config.muscle_path
        self.hmmer_path = run_config.hmmer_path

        self.msa_fasta = run_config.msa_fasta
        self.msa_name = run_config.msa_name
        self.hmm_model = run_config.hmm_model
        
    def MSA_create(self):
        if not os.path.exists(f'alignments/muscle_outputs//{self.msa_name}.afa'):
            print('Creating MSA using MUSCLE...')
            subprocess.run(f'{self.muscle_path} -align {self.msa_fasta} -output ../../alignments/muscle_outputs/{self.msa_name}.afa > ../../alignments/muscle_outputs/muscle_{self.msa_name}.log', shell=True, executable="/bin/zsh")

    def build(self):
        if not os.path.exists(f'hmm_model/{self.hmm_model}.hmm'):
            print('Building HMM model...')
            print(f'{self.hmmer_path}/hmmbuild --amino ../../hmm_model/{self.hmm_model}.hmm ../../alignments/muscle_outputs/{self.msa_name}.afa > ../../hmm_model/{self.hmm_model}_hmmbuild.log')
            subprocess.run(f'{self.hmmer_path}/hmmbuild --amino ../../hmm_model/{self.hmm_model}.hmm ../../alignments/muscle_outputs/{self.msa_name}.afa > ../../hmm_model/{self.hmm_model}_hmmbuild.log', shell=True, executable="/bin/zsh")

    def align(self,fasta_name,fasta_ext = 'faa',save_path = '../datasets'):
        subprocess.run(f'{self.hmmer_path}/hmmalign --trim --outformat afa ../../hmm_model/{self.hmm_model}.hmm ../datasets/{fasta_name}.{fasta_ext} > ../../alignments/{fasta_name}.afa', shell=True, executable="/bin/zsh")
        print(f'{self.hmmer_path}/hmmalign --trim --outformat afa ../../hmm_model/{self.hmm_model}.hmm ../datasets/{fasta_name}.{fasta_ext} > ../../alignments/{fasta_name}.afa')
        hmm_msa = AlignIO.read(f'../../alignments/{fasta_name}.afa', 'fasta')
        with open(f'{save_path}/{fasta_name}_fix.afa','w') as out_file:
            for prot in hmm_msa:
                out_file.write('>' + prot.id + '\n')
                out_file.write(helper_functions.fix_seqs(str(prot.seq)))
                out_file.write('\n')
        out_file.close()


class helper_functions:

    @staticmethod
    def fix_seqs(input_string):
        seq_ = input_string.replace('.', '')
        seq_ = ''.join([char for char in seq_ if not char.islower()])
        return seq_

    @staticmethod
    def gaplimiter(seq,threshold=0.7):
        gap_percent = seq.count('-') / len(seq)
        if gap_percent > threshold:
            return False
        else:
            return True
        
    @staticmethod
    def gaplimiter_df(df,threshold=0.7):
        return df[df['trimmed_afa'].apply(lambda x: gaplimiter(x,threshold))]

    @staticmethod
    def msa_column_extractor():
        with open(run.pickles+ "/msa_columns.pkl", "rb") as input_file:
            msa_columns = pickle.load(input_file)
        return msa_columns
    
    @staticmethod
    def apply_msa_mask(msa_df):
        msa_columns = msa_column_extractor()
        msa_df['trimmed_afa'] = msa_df['afa'].apply(lambda x: ''.join([x[i] for i in msa_columns]))
        return msa_df

    @staticmethod
    def aligned_df(afa_path):
        dict = {seq.description: fix_seqs(str(seq.seq)) for seq in SeqIO.parse(afa_path, "fasta")}
        df = pd.DataFrame.from_dict(dict, orient='index', columns=['afa'])
        df = apply_msa_mask(df)
        return df

    @staticmethod
    def encode_custom_seqs(afa_path,out_path,model= f'vae_fold_0.model',batch_size=1):
        df = aligned_df(afa_path)
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
    
    