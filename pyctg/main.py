import pandas as pd
import numpy as np
from synergy.utils.plots import plot_heatmap


class CTG_synergy:
    """Data class for CTG synergy analysis
    Assuming the experiment is a dose-titration experiment with 2 drugs
    happening in the middle 60 wells of a 96-well plate
    """
    def __init__(self, df, wide_treatment, narrow_treatment):
        self.df = df
        self.wide_treatment = wide_treatment
        self.narrow_treatment = narrow_treatment
    
    def ave_replicates(self):
        df = self.df.copy()
        df = df.set_index(['cell_type',self.narrow_treatment,self.wide_treatment]).pivot(
            columns='replicate', values='ctg'
        ).mean(axis=1).reset_index()
        df.columns = df.columns[:-1].tolist() + ['ctg']
        
        return df
    
    def extract_single_treatment(self, treatment_col):
        if treatment_col == self.wide_treatment:
            other_treatment_col = self.narrow_treatment
        elif treatment_col == self.narrow_treatment:
            other_treatment_col = self.wide_treatment
        else:
            raise ValueError(f"expected {self.wide_treatment} or {self.narrow_treatment}")

        df = self.df.copy()
        df = df.query(f'{other_treatment_col} == 0').drop(columns=[other_treatment_col])

        # Normalize each replicate's control (i.e. 0 treatment) for each cell type
        control_values = df[df[treatment_col] == 0].groupby(['cell_type', 'replicate'])['ctg'].mean().reset_index()
        control_values = control_values.rename(columns={'ctg': 'control_value'})
        df = df.merge(control_values, on=['cell_type', 'replicate'])
        df['ctg'] = df['ctg'] / df['control_value'] * 100
        df = df.drop(columns=['control_value'])

        df = df.rename(columns={treatment_col: 'Compound Conc'})
        df.insert(0, 'treatment', treatment_col)

        df.query('`Compound Conc` > 0', inplace=True)

        return df
    
    def fit_synergy_model(self):
        # model = MuSyC(E_bounds=(0, 1), h_bounds=(1e-3, 1e3), alpha_bounds=(1e-3, 1e3))
        # model = Bliss()
        # https://github.com/djwooten/synergy/issues/40
        # TODO: make sure how to normalize model fit for synergy values
        # Fit the model (bootstrap_iterations is an option for some models to estimate parameter
        # _ = model.fit(d1, d2, E)
        # S = model._get_synergy(d1, d2, E)
        pass

    def plot_synergy_heatmap(self, query, ax, xlabel='auto', ylabel='auto', remove_ticks=False, title=None):
        
        df = self.ave_replicates().query(query).copy()

        # Prepare the input data to be fit
        d1 = df[self.wide_treatment].to_numpy().astype(float)
        d2 = df[self.narrow_treatment].to_numpy().astype(float)
        E = df["ctg"].to_numpy().astype(float)

        if xlabel == 'auto':
            xlabel = self.wide_treatment
        if ylabel == 'auto':
            ylabel = self.narrow_treatment

        plot_heatmap(
            d1, d2, 
            E, 
            title=title,
            xlabel=f"\n{xlabel}", ylabel=f"{ylabel}\n",
            cmap="PRGn", 
            # center_on_zero=True,
            # vmin=-1, vmax=1,
            ax=ax,
        )
        if remove_ticks:
            # remove ticks and ticks bar
            ax.set_xticks([])
            ax.set_yticks([])
            ax.tick_params(axis='both', which='both', length=0)
            
            # remove color bar from plot
            ax.collections[0].colorbar.remove()

        # TODO: Use dose range from data to set x/y-ticks
        # ax.set_xticks(df.Idasanutlin.unique().round(decimals=2).astype(str).tolist())


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