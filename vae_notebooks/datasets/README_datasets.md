# Datasets for Variational autoencoder for Protein Ancestral Sequence Reconstruction

This directory contains all datasets described in the paper:<br>
Kohout P, Vasina M, Majerova M, Novakova V, Damborsky J, Bednar D, Marek M, Prokop Z, Mazurenko S. Engineering Dehalogenase Enzymes using Variational Autoencoder-Generated Latent Spaces and Microfluidics. (DOI 10.26434/chemrxiv-2023-jcds7)

- **PF00561_full.txt** : pfam dataset of AB hydrolases
- **hldI_IV.fa**  : HLDs dataset used for Model 1
- **hldI_II.fa**     : HLD I-II dataset used for Model 2-4

The other files are:

- **hldI_II_soluprot.txt** : predicted solubility data for HLD I-II dataset by EnzymeMiner
- **ancestors.fasta** : Babkova's ancestral sequences reconstructed via conventional strategies
- **Statistics/** : tree in Newick format by Fireprot-Asr, created for HDL I-II dataset only and already aligned. 
                    If you wish to measure tree mapping in latent space, create new trees by Fireprot Asr.
                    
How to prepare phylogenetic trees for validation of evolutionary latent space properties:

- **generate sampled MSA from training dataset for tree reconstruction** : you have to prepare MSA subsamples for phylogenetic tree reconstruction; there is a task scenario ready for that, 
    
    ```
    python3 runner.py run_task.py --run_stats_package_fireprot --json config_path.json
    ```
    this will prepare MSAs into *datasets/Statistics/data_source/fireprots_msas*
    
- **prepare trees**: use FireProtAsr <https://loschmidt.chemi.muni.cz/fireprotasr/> tool to prepare your trees. Use tab **User Data** and submit individual trees one by one.
- **update tree repository**: once Fireprot Jobs are ready, check that they contain ancestors (nothing failed) and place bigMSA**N**.fasta and tree.nwk (rename to msa_tree**N**.nwk) into *datasets/Statistics/tree_sources/<sub>your_dataset_name</sub>* directory. **N** refers to job sequential number (e.g. for job with MSA1.fasta you may name it as bigMSA1.fasta, msa_tree1.nwk)
- **run tree mapping task**:
    ```
    python3 runner.py run_task.py --run_stats_package_tree --json config_path.json
    ```
maps tree branches into the latent space and makes an evolutionary analysis (Figure 2 D-F in the paper)
