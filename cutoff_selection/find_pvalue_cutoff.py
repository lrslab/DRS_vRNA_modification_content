import pandas as pd
from plotnine import *
import numpy as np
from scipy import stats
import statistics
from matplotlib import pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
# After merge all output files from each method
path="SINV/subsets/"
differr_path=path+"differr.bed"
drummer_path=path+"drummer.bed"
eligos2_path=path+"eligos2.bed"
tombo_path=path+"tombo.bed"
nanocompore_path=path+"nanocompore.tsv"
xpore_path=path+"xpore.csv"
output_path="SINV/subsets/"
def plot_differr():
    df=pd.read_csv(differr_path,sep='\t',header=None)
    column_array = list(range(df.shape[1]))
    # df[15]=-np.log10(df[15].values)
    column_array[9]='value'
    df.columns=column_array
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)
    data =df['value'].values
    # identify_gaussian(data)

    data=sorted(data.tolist())
    cutoff=round(data[-len(data)//100])
    #绘制直方图并添加置信区间
    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        geom_density(color='black') +
        labs(title='Differr', x='-LogFDR', y='Frequency')+
        theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (4, 4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                title=element_text(size=13,
                                              family="Arial")
                        )
    )
    ggsave(filename=output_path+"differr.pdf",plot=histogram,dpi=300)
    print(histogram)

# g =  ggplot(df, aes(x='LogPvalue'))\
#      + geom_density()\
#      + labs(title='Distribution Plot(differr)', x='-LogPvalue', y='Frequency')\
#     + theme_bw()
# print(g)

def plot_drummer():
    df=pd.read_csv(drummer_path,sep='\t')
    # column_array = list(range(df.shape[1]))
    # # df[15]=-np.log10(df[15].values)
    # column_array[6]='value'
    df['value']= np.abs(np.log10(df['p_values_OR_adj']) * (-1))
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)
    data =df['value'].values
    # identify_gaussian(data)

    data=sorted(data.tolist())
    cutoff=round(data[-len(data)//100])
    print(cutoff)

    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        geom_density(color='black') +
        labs(title='DRUMMER', x='-LogPvalue', y='Frequency')+
        theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (4, 4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                title=element_text(size=13,
                                   family="Arial")
                        )
    )
    ggsave(filename=output_path+"drummer.pdf",plot=histogram,dpi=300)
    print(histogram)

def plot_eligos2():
    df=pd.read_csv(eligos2_path,sep='\t',header=None,skiprows=1)
    column_array = list(range(df.shape[1]))
    df[16]=-np.log10(df[16].values)
    column_array[16]='value'
    df.columns=column_array
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)
    data =df['value'].values
    # identify_gaussian(data)

    data=sorted(data.tolist())
    cutoff=round(data[-len(data)//100])
    print(cutoff)
    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        geom_density(color='black') +
        labs(title='ELIGOS2', x='-LogPvalue', y='Frequency')+
        geom_vline(xintercept=cutoff, linetype='dashed', color='red') +
        theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (4, 4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                title=element_text(size=13,
                                   family="Arial")
                        )
    )
    ggsave(filename=output_path+"eligos2.pdf",plot=histogram,dpi=300)
    print(histogram)

def plot_tombo():
    df=pd.read_csv(tombo_path,sep='\t',header=None)
    column_array = list(range(df.shape[1]))
    # df[15]=np.log2(df[15].values)
    column_array[6]='value'
    df.columns=column_array
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)
    data =df['value'].values
    # identify_gaussian(data)

    data=sorted(data.tolist())
    cutoff=round(data[-len(data)//100])
    print(cutoff)
    #绘制直方图并添加置信区间
    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        geom_density(color='black') +
        labs(title='Tombo', x='-LogPvalue', y='Frequency')+
        geom_vline(xintercept=cutoff, linetype='dashed', color='red') +
        theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (4, 4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                title=element_text(size=13,
                                   family="Arial")
                        )
    )
    ggsave(filename=output_path+"tombo.pdf",plot=histogram,dpi=300)
    print(histogram)


def plot_nanocompore():
    df=pd.read_csv(nanocompore_path,sep='\t',header=None,skiprows=1)
    column_array = list(range(df.shape[1]))
    df[6]=-np.log10(df[6].values)
    column_array[6]='value'
    df.columns=column_array
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)

    data = df['value'].values
    data=sorted(data.tolist())
    cutoff=round(data[-len(data)//100])
    print(cutoff)
    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        geom_density(color='black') +
        labs(title='Nanocompore', x='-LogPvalue', y='Frequency')+
        theme_bw()\
        +xlim(0,1)\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (4, 4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                title=element_text(size=13,
                                   family="Arial")
                        )
    )
    ggsave(filename=output_path+"nanocompore.pdf",plot=histogram,dpi=300)
    print(histogram)

def plot_xpore():
    df=pd.read_csv(xpore_path,header=None,skiprows=1)
    column_array = list(range(df.shape[1]))
    df[4]=-np.log10(df[4].values)
    column_array[4]='value'
    df.columns=column_array
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)
    data =df['value'].values


    data=sorted(data.tolist())
    cutoff=round(data[-len(data)//100])
    print(cutoff)
    #绘制直方图并添加置信区间
    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        geom_density(color='black') +
        geom_vline(xintercept=cutoff, linetype='dashed', color='red') +
        labs(title='Xpore', x='-LogPvalue', y='Frequency')+
        theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (4, 4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                title=element_text(size=13,
                                   family="Arial")
                        )
    )
    ggsave(filename=output_path+"xpore.pdf",plot=histogram,dpi=300)
    print(histogram)


# plot_tombo()
#
# plot_xpore()
# plot_drummer()
# plot_differr()
# plot_eligos2()
plot_nanocompore()
