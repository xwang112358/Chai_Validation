from pathlib import Path
from pdb2fasta import extract_chain_A_fasta
import numpy as np
import torch
import os
from chai_lab.chai1 import run_inference
import argparse


# def sort_af2pdb_list(pdb_list):
#     return sorted(pdb_list, key=lambda x: (int(x.split('_')[3]), int(x.split('_')[-])))

# We use fasta-like format for inputs.
# - each entity encodes protein, ligand, RNA or DNA
# - each entity is labeled with unique name;
# - ligands are encoded with SMILES; modified residues encoded like AAA(SEP)AAA

# Example given below, just modify it

# add argument parser

cd20_fasta = """
>protein|cd20-chain-A
MRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAG
>protein|cd20-chain-B
MRESKTLGAVQIMNGLFHIALGGLLMIPAGIYAPICVTVWYPLWGGIMYIISGSLLAATEKNSRKCLVKGKMIMNSLSLFAAISGMILSIMDILNIKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQYCYSIQSLFLGILSVMLIFAFFQELVIAG
""".strip()

parser = argparse.ArgumentParser()
parser.add_argument("--name", type=str, required=True, help="name of the design")
args = parser.parse_args()
name = args.name

pdb_dir = f"./final_filtered/{name}"
pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]

# sort pdb files by string
pdb_files = sorted(pdb_files)

results = {
    'aggregate_score': [],
    'ptm': [],
    'iptm': [],
    'per_chain_ptm': [],
    'per_chain_pair_iptm': [],
    'has_inter_chain_clashes': [],
    'chain_chain_clashes': [],
    'pdb_name': []
}

for pdb_file in pdb_files:
    pdb_name = pdb_file.split(".")[0]

    if os.path.exists(f"./outputs/{pdb_name}"):
        print(f"Skipping {pdb_file} because it already exists in outputs")
        continue

    fasta_path = Path(f"/home/allenwang/protein/chai-lab/inputs/{pdb_name}.fasta")
    print(pdb_file)
    binder = extract_chain_A_fasta(os.path.join(pdb_dir, pdb_file)).strip()
    complex_fasta = cd20_fasta + "\n" + binder

    fasta_path.write_text(complex_fasta)

    output_dir = Path(f"/home/allenwang/protein/chai-lab/outputs/{pdb_name}") 

    candidates = run_inference(
        fasta_file=fasta_path,
        output_dir=output_dir,
        # 'default' setup
        num_trunk_recycles=3,
        num_diffn_timesteps=200,
        seed=42,
        device=torch.device("cuda:0"),
        use_esm_embeddings=True,
    )
    
    scores = np.load(output_dir.joinpath("scores.model_idx_0.npz"))
    results['aggregate_score'].append(scores['aggregate_score'])
    results['ptm'].append(scores['ptm'])
    results['iptm'].append(scores['iptm'])
    results['per_chain_ptm'].append(scores['per_chain_ptm'])
    results['per_chain_pair_iptm'].append(scores['per_chain_pair_iptm'])
    results['has_inter_chain_clashes'].append(scores['has_inter_chain_clashes'])
    results['chain_chain_clashes'].append(scores['chain_chain_clashes'])
    results['pdb_name'].append(pdb_file)


# save results as npz
np.savez(f'results/{name}_results_2.npz', **results)


# cif_paths = candidates.cif_paths
# scores = [rd.aggregate_score for rd in candidates.ranking_data]


# # Load pTM, ipTM, pLDDTs and clash scores for sample 2
# scores = np.load(output_dir.joinpath("scores.model_idx_2.npz"))
