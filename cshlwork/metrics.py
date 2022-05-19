#!/usr/bin/env python
#
# Various statistics, plotting utils
#
#
import logging
import pandas as pd
import numpy as np
import seaborn as sns
import math
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px


def calc_auroc(df, p_col, t_col ):
    t_array = np.array(df.sort_values( by=p_col, ascending=False).val)
    t_val = t_array.cumsum() / t_array.sum()
    slicewidth = 1 / (len(t_val) ) 
    auroc = 0.0
    
    for i,v in enumerate(t_val):
        auroc = auroc + (v * slicewidth)
        try:
            triarea = ( slicewidth * abs(t_val[i] - t_val[i+1]) ) / 2 
            auroc = auroc +  triarea
        except IndexError:
            pass
    auroc_str = "{0:.3g}".format(auroc)
    return auroc_str


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return zip(*[iter(iterable)]*n)

def plot_auroc(df, pcol, vcol, diagonal=False, title=None):
    '''
    col_list should be a list of columns in df, where first is model probabilities, and 
    second is one-hot membership in class (1=member,0=non-member). 
    
    '''
    auroc_str = "NA"
    #for pcol, vcol in grouped(col_list, 2):
    logging.debug(f"prob_col={pcol}, valid_col={vcol}")

    gdf = df.sort_values( by=pcol, ascending=False)       
    gdf['fval'] = 1 - df[vcol]
    gdf['t_val']  = gdf[vcol].cumsum() / gdf[vcol].sum()
    gdf['f_val'] = gdf['fval'].cumsum() / gdf['fval'].sum()
    num_items = len(gdf['t_val'])
    
    auroc_str = calc_auroc(gdf, pcol, vcol)
            
    palette = sns.color_palette("mako_r", 6)
    ax = sns.lineplot(data=gdf, 
                      x='f_val', 
                      y='t_val',
                     palette = palette,
                     #legend = 'full',
                     label=f'AUROC={auroc_str} n={num_items}'
                     )
    #ax.text(0.65,0.65,f'AUROC={auroc_str}', fontsize=14)
    #ax.add_shape(type='line', line=dict(dash='dash'), x0=0, x1=1, y0=0, y1=1)
    if title is not None:
        ax.set(title=title)
    else:    
        ax.set(title=f'AUROC {pcol} vs. {vcol}')
    ax.set(ylabel='TPR',xlabel='FPR')
    return ax


def plot_auroc_group(df, pcol, vcol, diagonal=False, title=None, group_col=None):
    '''
    col_list should be a list of columns in df, where first is model probabilities, and 
    second is one-hot membership in class (1=member,0=non-member). 
    
    '''
    auroc_str = "NA"
           
    logging.debug(f"prob_col={pcol}, valid_col={vcol} group_col={group_col}")
    #for grp in list(df[group_col].unique() ):
    #    gdf = df[df[group_col] == grp]
    #    logging.debug(f'group df len={len(gdf)}')
    
    gdf = df.sort_values( by=pcol, ascending=False)       
    gdf['fval'] = 1 - df[vcol]
    gdf['t_val']  = gdf[vcol].cumsum() / gdf[vcol].sum()
    gdf['f_val'] = gdf['fval'].cumsum() / gdf['fval'].sum()
    num_items = len(gdf['t_val'])
    
    auroc_str = calc_auroc(gdf, pcol, vcol)
            
    #palette = sns.color_palette("mako_r", 6)
    ax = sns.lineplot(data=gdf, 
                      x='f_val', 
                      y='t_val',
                      hue=group_col,
                      #palette = palette,
                      legend = 'full',
                      #label=f'{pcol} AUROC = {auroc_str} n={num_items}'
                     )
    #ax.text(0.65,0.65,f'AUROC={auroc_str}', fontsize=14)
    #ax.add_shape(type='line', line=dict(dash='dash'), x0=0, x1=1, y0=0, y1=1)
    if title is not None:
        ax.set(title=title)
    else:    
        ax.set(title=f'AUROC {pcol} vs. {vcol}')
    ax.set(ylabel='TPR',xlabel='FPR')
    return ax

    

def plot_interactive_auroc(df, pcol, vcol):
    #threshold column is pcol
    tcol = pcol
    fig = px.area(
        data_frame=df, 
        x=fpr, 
        y=tpr,
        hover_data=thresholds, 
        title='ROC Curve'
    )
    fig.update_layout(
        autosize=False,
        width=500,
        height=500,
        margin=dict(l=30, r=30, b=30, t=30, pad=4),
        title_x=.5, # Centre title
        hovermode = 'closest',
        xaxis=dict(hoverformat='.4f'),
        yaxis=dict(hoverformat='.4f')
    )
    hovertemplate = 'False Positive Rate=%{x}<br>True Positive Rate=%{y}<br>Threshold=%{customdata[0]:.4f}<extra></extra>'
    fig.update_traces(hovertemplate=hovertemplate)
    
    # Add dashed line with a slope of 1
    fig.add_shape(type='line', line=dict(dash='dash'), x0=0, x1=1, y0=0, y1=1)
    fig.show()
    
    #                       plot_interactive_roc_curve(df=inputs, 
    #                       fpr='false_positive_rate', 
    #                       tpr='true_positive_rate', 
    #                       thresholds=['threshold'])
    