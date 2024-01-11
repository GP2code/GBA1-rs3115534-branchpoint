import pandas as pd

# Load the CSV file into a DataFrame
df = pd.read_csv('HBCC_82040_types.csv', header=None)
df.columns = ['Sequence', 'Cell_Type']

# Extract the sequence without the '-1' at the end
df['Sequence'] = df['Sequence'].str.rstrip('-1')

# Iterate through unique cell types
unique_cell_types = df['Cell_Type'].unique()
for cell_type in unique_cell_types:
    # Create a DataFrame for the current cell type
    cell_type_df = df[df['Cell_Type'] == cell_type]
    
    # Extract only the 'Sequence' column and save it to a text file
    output_file = f'{cell_type}.txt'
    cell_type_df['Sequence'].to_csv(output_file, header=False, index=False, sep='\t')
    