import pandas as pd

def load_data(drugs_data_path, targets_data_path):
        """
        Load drug and target data from CSV files and create lookup dictionaries.
        
        Creates two dictionaries:
        1. drug_mapping_dict: maps each drug to its set of protein targets
        2. fingerprint_dict: maps each drug to its molecular fingerprint set
        
        Inputs:
            None (uses paths stored in instance variables)
        
        Returns:
            None (stores data in instance variables)
        """
        # Load CSV files into pandas DataFrames
        drugs_df = pd.read_csv(drugs_data_path)
        targets_df = pd.read_csv(targets_data_path)
        
        # Create dictionary mapping each drug to set of protein targets it binds
        # Key: db_id (drug), Value: set of uniprot_accession (proteins)
        drug_mapping_dict = {}
        for _, row in targets_df.iterrows():
            drug_id = row['db_id']
            protein = row['uniprot_accession']
            
            if drug_id not in drug_mapping_dict:
                drug_mapping_dict[drug_id] = set()
            drug_mapping_dict[drug_id].add(protein)
        
        # Create dictionary mapping each drug to its fingerprint set
        # Key: db_id (drug), Value: set of fingerprint features (integers)
        fingerprint_dict = {}
        for _, row in drugs_df.iterrows():
            drug_id = row['db_id']
            # Parse space-separated fingerprint string into set of integers
            fingerprint = set(map(int, row['maccs'].split()))
            fingerprint_dict[drug_id] = fingerprint
        

        return drugs_df, targets_df, drug_mapping_dict, fingerprint_dict


def tanimoto_coeff(drug_a, drug_b, fingerprint_dict):
    """
    Calculate Tanimoto coefficient between two drugs based on fingerprints.
    
    The Tanimoto coefficient (Jaccard Index) measures chemical similarity:
    Tc = |fpt(A) ∩ fpt(B)| / |fpt(A) ∪ fpt(B)|
    
    Where fpt(A) and fpt(B) are the fingerprint feature sets.
    
    Inputs:
        drug_a: string, db_id of first drug
        drug_b: string, db_id of second drug
    
    Returns:
        float, Tanimoto coefficient between 0.0 (no similarity) and 1.0 (identical)
    """
    # Get fingerprint sets for both drugs
    fpt_a = fingerprint_dict[drug_a]
    fpt_b = fingerprint_dict[drug_b]
    
    # Calculate size of intersection (common features)
    intersection = len(fpt_a & fpt_b)
    # Calculate size of union (all unique features)
    union = len(fpt_a | fpt_b)
    
    # Return ratio, handling edge case of empty union
    return intersection / union if union > 0 else 0.0