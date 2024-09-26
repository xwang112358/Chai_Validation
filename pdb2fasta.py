from Bio import PDB
from Bio.SeqUtils import seq1

def extract_chain_A_fasta(pdb_file_path):
    """
    Extracts the sequence of chain A from a PDB file and returns it as a FASTA string.

    Args:
    - pdb_file_path (str): Path to the PDB file.

    Returns:
    - str: A FASTA-formatted string with the chain A sequence.
    """
    # Initialize the PDB parser
    parser = PDB.PDBParser(QUIET=True)
    
    # Parse the PDB file
    structure = parser.get_structure('protein', pdb_file_path)
    
    # Extract the sequence of chain A
    chain_A = None
    for model in structure:
        for chain in model:
            if chain.id == 'A':  # Select chain A
                chain_A = chain
                break
        if chain_A:
            break
    
    # Ensure chain A was found
    if chain_A is None:
        raise ValueError("Chain A not found in the PDB file.")
    
    # Extract the sequence in one-letter code
    sequence = ""
    for residue in chain_A:
        if PDB.is_aa(residue, standard=True):  # Only standard amino acids
            sequence += seq1(residue.resname)
    
    # Format the sequence as a FASTA string
    fasta_string = f">protein|binder\n{sequence}"
    
    return fasta_string

# # Example usage
# pdb_file_path = 'example.pdb'  # Replace with your actual PDB file path
# fasta_string = extract_chain_A_fasta(pdb_file_path)
# print(fasta_string)
