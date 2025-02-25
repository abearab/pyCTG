import pandas as pd


def read_CTG_titration_data(filename):
    data = pd.read_csv(
        filename,
        sep='\t',
        header=0, 
        # index_col=[0, 1]
    )
    
    df = data.melt(
        id_vars=['treatment','cell_type','replicate'], 
        value_name='ctg', 
        var_name='Compound Conc'
    )

    df['Compound Conc'] = df['Compound Conc'].astype(float) * 1000 # convert to nM
    
    # round to 3 decimals
    for col in ['ctg', 'Compound Conc']:
        df[col] = df[col].round(decimals=3)
    
    return df
