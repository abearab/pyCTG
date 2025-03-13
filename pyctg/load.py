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


def read_CTG_synergy_data(filename):
    data = pd.read_csv(filename,sep='\t', header=0, index_col=None, skiprows=1)

    #TODO: come up with a better way to get the treatment names
    wide_treatment = pd.read_csv(filename,sep='\t', nrows=1, header=None).T.dropna().values[0][0]
    narrow_treatment = data.columns[2]

    df = data.melt(
        id_vars=['cell_type','replicate',narrow_treatment], 
        value_name='ctg', 
        var_name=wide_treatment
    )
    df[wide_treatment] = df[wide_treatment].astype(float) * 1000 # convert to nM
    df[narrow_treatment] = df[narrow_treatment].astype(float) * 1000 # convert to nM

    # round to 3 decimals
    for col in ['ctg', narrow_treatment, wide_treatment]:
        df[col] = df[col].round(decimals=3)

    return CTG_synergy(df, wide_treatment, narrow_treatment)