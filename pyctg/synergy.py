import pandas as pd
import numpy as np
from synergy.combination import Bliss
from synergy.combination.loewe import Loewe
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
    
    def extract_single_treatment(self, treatment_col, treatment_dose=0):
        if treatment_col == self.wide_treatment:
            other_treatment_col = self.narrow_treatment
        elif treatment_col == self.narrow_treatment:
            other_treatment_col = self.wide_treatment
        else:
            raise ValueError(f"expected {self.wide_treatment} or {self.narrow_treatment}")

        df = self.df.copy()
        df = df.query(f'{other_treatment_col} == {treatment_dose}').drop(columns=[other_treatment_col])

        df = df.rename(columns={treatment_col: 'Compound Conc'})
        df.insert(0, 'treatment', treatment_col)

        df.query('`Compound Conc` > 0', inplace=True)

        return df
    
    def plot_synergy_heatmap(self, query, ax, value_col='viability', xlabel='auto', ylabel='auto', remove_ticks=False, title=None, cmap="PRGn", colorbar=True, **args):
        
        # calculate bliss synergy if needed
        if value_col in ['bliss','loewe'] and value_col not in self.df.columns:
            df = self.calculate_synergy(method=value_col, inplace=False)

        df = self._ave_replicates(value_col=value_col).query(query).copy()

        # Prepare the input data to be fit
        d1 = df[self.wide_treatment].to_numpy().astype(float)
        d2 = df[self.narrow_treatment].to_numpy().astype(float)
        E = df[value_col].to_numpy().astype(float)

        if xlabel == 'auto':
            xlabel = self.wide_treatment
        if ylabel == 'auto':
            ylabel = self.narrow_treatment

        plot_heatmap(
            d1, d2, 
            E, 
            title=title,
            xlabel=f"\n{xlabel}", ylabel=f"{ylabel}\n",
            cmap=cmap,
            # center_on_zero=True,
            # vmin=-1, vmax=1,
            ax=ax,
            **args
        )
        if remove_ticks:
            # remove ticks and ticks bar
            ax.set_xticks([])
            ax.set_yticks([])
            ax.tick_params(axis='both', which='both', length=0)
        
        if not colorbar:
            # remove color bar from plot
            ax.collections[0].colorbar.remove()

        # TODO: Use dose range from data to set x/y-ticks
        # ax.set_xticks(df.Idasanutlin.unique().round(decimals=2).astype(str).tolist())

    def calculate_synergy(self, method='bliss', inplace=True, **kwargs):
        #TODO add checks here
        # if method not in ['bliss', 'loewe']:

        df = self.df.copy()

        # https://github.com/djwooten/synergy/issues/40
        # TODO: make sure how to normalize model fit for synergy values
        
        df[method] = np.nan

        for _,row in df.query(
            f'`{self.wide_treatment}` == 0 & `{self.narrow_treatment}` == 0').iterrows():
            single_plate = df.loc[
                (df.cell_type == row['cell_type']) & 
                (df.replicate == row['replicate']),:]

            if method == 'bliss':
                model = Bliss()
            elif method == 'loewe':
                model = Loewe(mode="delta_hsa")
            
            res = model.fit(
                d1=single_plate[self.wide_treatment].to_numpy(), 
                d2=single_plate[self.narrow_treatment].to_numpy(), 
                E=single_plate['viability'].to_numpy(),
                **kwargs
            )
            
            df.loc[single_plate.index, method] = res
        
        if inplace:
            self.df = df
        else:
            return df        
    
    def _ave_replicates(self, value_col):
        df = self.df.copy()
        df = df.set_index(['cell_type',self.narrow_treatment,self.wide_treatment]).pivot(
            columns='replicate', values=value_col
        ).mean(axis=1).reset_index()
        df.columns = df.columns[:-1].tolist() + [value_col]
        
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

    # calculate relative CTG values (normalized to baseline, i.e. no treatment)
    df['baseline'] = np.nan

    for _,row in df.query(
        f'`{wide_treatment}` == 0 & `{narrow_treatment}` == 0').iterrows():
        df.loc[
            (df.cell_type == row['cell_type']) & 
            (df.replicate == row['replicate']), 
            'baseline'] = row['ctg']

    # viability
    df['viability'] = (df['ctg'] / df['baseline'])
    del df['baseline']

    return CTG_synergy(df, wide_treatment, narrow_treatment)
