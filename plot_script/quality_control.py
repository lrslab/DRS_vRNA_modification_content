import pandas as pd
import numpy as np
from plotnine import *
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
def extract_mis_ratio(df,mis,label):
    L_rep1 = pd.read_csv(df, sep='\t')
    if L_rep1.shape[0]>1000000:
        L_rep1=L_rep1.sample(frac=0.2)
    print(L_rep1.shape)
    # L_rep1=L_rep1[L_rep1['del']<1000]
    # L_rep1=L_rep1[L_rep1['ins']<500]
    print(L_rep1.shape)
    if mis=='acc' or mis == 'iden' or mis == 'mapq':
        L_rep1 = L_rep1[[mis]]
        L_rep1["sample"] = label
        df = L_rep1.sort_values(by=mis)
        # 删除最xiao的20个数
        top_num = int(df.shape[0] * 0.01)
        df = df.drop(df.index[:top_num].values, axis=0)
        return df
    elif mis=='read_length':
        L_rep1 = L_rep1[["read_length"]]
        L_rep1["sample"] = label
        df = L_rep1.sort_values(by=mis)
        # 删除最大的20个数
        top_num = int(df.shape[0] * 0.005)
        df = df.drop(df.index[:top_num].values, axis=0)
        df = df.drop(df.index[-top_num:].values, axis=0)
        print(df.shape)

    elif mis=='ratio':
        temp=(L_rep1['aligned_ref_len'].sum(), L_rep1['length'].sum())
        return temp

    elif mis=='ins' or mis=='sub' or mis=='del' :
        temp=(L_rep1[mis]/ L_rep1['length'])
        L_rep1[mis] = temp
        L_rep1 = L_rep1[[mis]]
        L_rep1["sample"] = label

        df = L_rep1.sort_values(by=mis)
        # 删除最大的20个数
        top_num=int(df.shape[0]*0.01)
        df = df.drop(df.index[-top_num:].values, axis=0)
        return  df

    return df

def draw_virus():

    wt_path = "SINV/BHK_polyA/subsample_gRNA/stats_filtered.txt"
    ivt_path = "C636_polyA/subsample_gRNA/stats_filtered.txt"
    ivt_path_2 = "IVT/subsample_gRNA/stats_filtered.txt"
    col = 'mapq'
    title={'sub':"Mismatch",
           'ins':"Insertion",
           'del':'Deletion',
           'iden':"Identity",
           'acc':"Accuracy",
           "read_length":"Read length",
           'mapq':'mapq'}
    wt_path = extract_mis_ratio(wt_path, col, 'BHK')
    ivt = extract_mis_ratio(ivt_path, col, 'C636')
    results_df = pd.concat([wt_path, ivt], sort=False)
    ivt = extract_mis_ratio(ivt_path_2, col, 'IVT')
    results_df = pd.concat([results_df, ivt], sort=False)
    category = pd.api.types.CategoricalDtype(categories=["BHK","C636","IVT"], ordered=True)
    results_df['sample'] = results_df['sample'].astype(category)

    tmp_plot = ggplot(results_df, aes(x='sample', y=col,fill='sample')) + \
                geom_violin(color='none') + \
                scale_fill_manual(values=["#8ECFC9",'#FFBE7A', '#FA7F6F']) + \
                geom_boxplot(width=0.1,outlier_alpha=0) + \
                theme_bw() \
               + theme(legend_key_height=10,
                       legend_key_width=5,
                       figure_size=(4, 4),
                       panel_grid_minor=element_blank(),
                       axis_text=element_text(size=13,family='Arial'),
                       axis_title=element_text(size=13,family='Arial'),
                       legend_position='none'
                       ) \
            +labs( x='Sample', y=title[col])
    ggsave(filename="quality_control/"+col+".pdf",plot=tmp_plot,dpi=300)
    print(tmp_plot)
    print(1)