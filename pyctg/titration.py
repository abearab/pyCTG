import pandas as pd
import numpy as np
from py50 import Calculator, PlotCurve, CBMARKERS, CBPALETTE

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
    
    #TODO: calculate relative CTG values (normalized to baseline, i.e. no treatment) 

    return df


def plot_CTG_titration(ctg_data, treatment_name, value_col='viability', title=None, ymax=115, ymin=-10):
    df = ctg_data.query(
        f'treatment == "{treatment_name}"' # & `Compound Conc` < 10000'
    ).drop(columns=['treatment']).pivot(index=['cell_type','Compound Conc'], columns='replicate', values=value_col)
    df.reset_index(inplace=True)
    df.columns.name = None

    data = Calculator(df)

    ic50 = data.calculate_ic50(
        name_col='cell_type',
        concentration_col='Compound Conc',
        response_col=['rep1', 'rep2', 'rep3'],
    )

    data = PlotCurve(df)
    
    treatment_label = f"[{treatment_name}] (nM)"
    if '+' in treatment_name:
        treatment_label = treatment_label.replace('+', '],[')

    figure = data.multi_curve_plot(name_col='cell_type',
                                concentration_col='Compound Conc',
                                response_col=['rep1', 'rep2', 'rep3'], # % Inhibition Avg
                                plot_title=title,
                                xlabel=treatment_label,
                                ylabel='Viability %',
                                errorbar='sd',
                                legend=True,
                                #    hline=75,
                                ymax=ymax, ymin=ymin,
                                line_color=CBPALETTE, marker=CBMARKERS
    )

    return figure, ic50
