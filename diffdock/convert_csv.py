import pandas as pd

input_df = pd.read_csv('SMILES_inh.csv')

# input_df.sort_values(by='QED_Filter', ascending=False)

input_df = input_df[['inhibitor_id', 'smiles']]
input_df = input_df.iloc[:1000]

input_df.reset_index(inplace=True)
input_df.rename(columns={'inhibitor_id': 'complex_name', 'smiles': 'ligand_description'}, inplace=True)
input_df['protein_path'] = 'data/ipz/4ypx.pdb'
input_df['protein_sequence'] = None
input_df = input_df[['complex_name','protein_path','ligand_description','protein_sequence']]
input_df.to_csv('1000_samples.csv', index=None)
