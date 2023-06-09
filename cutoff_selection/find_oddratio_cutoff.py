import pandas as pd
from plotnine import *
import numpy as np
from scipy import stats
import statistics

from matplotlib import pyplot as plt
# After merge all output files from each method

plt.rcParams['pdf.fonttype'] = 42

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
    column_array[6]='value'
    df.columns=column_array
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)
    data =df['value'].values
    # identify_gaussian(data)


    confidence_level = 0.98
    sample_mean = np.mean(data)  # 使用样本数据计算平均值
    standard_error = np.std(data)  # 使用样本数据计算标准误差
    degrees_freedom = len(data) - 1  # 计算自由度，这里用 n-1

    # 使用 t 分布计算置信区间，对应的函数是 t.interval()
    lower_bound, upper_bound = stats.norm.interval(confidence_level, loc=sample_mean, scale=standard_error)

    sorted_data = sorted(data.tolist())
    midpoint = sorted_data[len(sorted_data)-(len(sorted_data) // 200)]

    # 计算一百分位数

    print(upper_bound)
    #绘制直方图并添加置信区间
    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        # geom_vline(xintercept=lower_bound, linetype='dashed', color='red') +
        geom_vline(xintercept=upper_bound, linetype='dashed', color='red') +
        geom_density(color='black') +
        labs(title='Differr', x='LogOddRatio', y='Frequency')+
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
    ggsave(filename=output_path+"oddratio/differr_3.31.pdf",plot=histogram,dpi=300)

def plot_drummer():
    df=pd.read_csv(drummer_path,sep='\t')
    # column_array = list(range(df.shape[1]))
    # # df[15]=-np.log10(df[15].values)
    # column_array[6]='value'
    df['value']= np.log2(df["odds_ratio"].values)
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)
    data =df['value'].values
    # identify_gaussian(data)


    confidence_level = 0.98
    sample_mean = np.mean(data)  # 使用样本数据计算平均值
    standard_error = np.std(data)  # 使用样本数据计算标准误差
    degrees_freedom = len(data) - 1  # 计算自由度，这里用 n-1

    # 使用 t 分布计算置信区间，对应的函数是 t.interval()
    lower_bound, upper_bound = stats.norm.interval(confidence_level, loc=sample_mean, scale=standard_error)

    sorted_data = sorted(data.tolist())
    midpoint = sorted_data[len(sorted_data)-(len(sorted_data) // 200)]

    # 计算一百分位数

    print(upper_bound)
    #绘制直方图并添加置信区间
    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        # geom_vline(xintercept=lower_bound, linetype='dashed', color='red') +
        geom_vline(xintercept=upper_bound, linetype='dashed', color='red') +
        geom_density(color='black') +
        labs(title='DRUMMER', x='LogOddRatio', y='Frequency')+
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
    ggsave(filename=output_path+"oddratio/drummer_1.10.pdf",plot=histogram,dpi=300)
# g =  ggplot(df, aes(x='LogPvalue'))\
#      + geom_density()\
#      + labs(title='Distribution Plot(differr)', x='-LogPvalue', y='Frequency')\
#     + theme_bw()
# print(g)

def plot_eligos2():
    df=pd.read_csv(eligos2_path,sep='\t',header=None,skiprows=1)
    column_array = list(range(df.shape[1]))
    df[15]=np.log2(df[15].values)
    column_array[15]='value'
    df.columns=column_array
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)
    data =df['value'].values


    confidence_level = 0.98
    sample_mean = np.mean(data)  # 使用样本数据计算平均值
    standard_error = np.std(data)  # 使用样本数据计算标准误差
    degrees_freedom = len(data) - 1  # 计算自由度，这里用 n-1

    # 使用 t 分布计算置信区间，对应的函数是 t.interval()
    lower_bound, upper_bound = stats.norm.interval(confidence_level, loc=sample_mean, scale=standard_error)

    sorted_data = sorted(data.tolist())
    midpoint = sorted_data[len(sorted_data)-(len(sorted_data) // 200)]

    # 计算一百分位数

    print(upper_bound)
    #绘制直方图并添加置信区间
    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        # geom_vline(xintercept=lower_bound, linetype='dashed', color='red') +
        geom_vline(xintercept=upper_bound, linetype='dashed', color='red') +
        geom_density(color='black') +
        labs(title='ELIGOS2', x='LogOddRatio', y='Frequency')+
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
    ggsave(filename=output_path+"oddratio/eligos2_1.01.pdf",plot=histogram,dpi=300)

def plot_tombo():
    df=pd.read_csv(tombo_path,sep='\t',header=None)
    column_array = list(range(df.shape[1]))
    # df[15]=np.log2(df[15].values)
    column_array[7]='value'
    df.columns=column_array
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)
    data =df['value'].values


    confidence_level = 0.99
    sample_mean = np.mean(data)  # 使用样本数据计算平均值
    standard_error = np.std(data)  # 使用样本数据计算标准误差
    degrees_freedom = len(data) - 1  # 计算自由度，这里用 n-1

    # 使用 t 分布计算置信区间，对应的函数是 t.interval()
    lower_bound, upper_bound = stats.norm.interval(confidence_level, loc=sample_mean, scale=standard_error)

    sorted_data = sorted(data.tolist())
    midpoint = sorted_data[len(sorted_data)-(len(sorted_data) // 200)]

    # 计算一百分位数

    print(upper_bound)
    #绘制直方图并添加置信区间
    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        geom_vline(xintercept=lower_bound, linetype='dashed', color='red') +
        geom_vline(xintercept=upper_bound, linetype='dashed', color='red') +
        geom_density(color='black') +
        labs(title='Tombo', x='Difference', y='Frequency')+
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
    ggsave(filename=output_path+"oddratio/tombo_0.06.pdf",plot=histogram,dpi=300)


def plot_nanocompore():
    df=pd.read_csv(nanocompore_path,sep='\t',header=None,skiprows=1)
    column_array = list(range(df.shape[1]))
    # df[15]=np.log2(df[15].values)
    column_array[-1]='value'
    df.columns=column_array
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)
    data =df['value'].values


    confidence_level = 0.99
    sample_mean = np.mean(data)  # 使用样本数据计算平均值
    standard_error = np.std(data)  # 使用样本数据计算标准误差
    degrees_freedom = len(data) - 1  # 计算自由度，这里用 n-1

    # 使用 t 分布计算置信区间，对应的函数是 t.interval()
    lower_bound, upper_bound = stats.norm.interval(confidence_level, loc=sample_mean, scale=standard_error)

    sorted_data = sorted(data.tolist())
    midpoint = sorted_data[len(sorted_data)-(len(sorted_data) // 200)]

    # 计算一百分位数

    print(upper_bound)
    #绘制直方图并添加置信区间
    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        geom_vline(xintercept=lower_bound, linetype='dashed', color='red') +
        geom_vline(xintercept=upper_bound, linetype='dashed', color='red') +
        geom_density(color='black') +
        xlim(-3, 3) +
        labs(title='Nanocompore', x='LogOddRatio', y='Frequency')+
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
    ggsave(filename=output_path+"oddratio/nanocompore_0.5.pdf",plot=histogram,dpi=300)

def plot_xpore():
    df=pd.read_csv(xpore_path,header=None,skiprows=1)
    column_array = list(range(df.shape[1]))
    # df[15]=np.log2(df[15].values)
    column_array[3]='value'
    df.columns=column_array
    df=df[['value']]
    df=df.replace([np.inf, -np.inf], np.nan).dropna()
    df=df.replace("NC", np.nan).dropna()
    df['value']=df['value'].astype(float)
    data =df['value'].values


    confidence_level = 0.99
    sample_mean = np.mean(data)  # 使用样本数据计算平均值
    standard_error = np.std(data)  # 使用样本数据计算标准误差
    degrees_freedom = len(data) - 1  # 计算自由度，这里用 n-1

    # 使用 t 分布计算置信区间，对应的函数是 t.interval()
    lower_bound, upper_bound = stats.norm.interval(confidence_level, loc=sample_mean, scale=standard_error)

    sorted_data = sorted(data.tolist())
    midpoint = sorted_data[len(sorted_data)-(len(sorted_data) // 200)]

    # 计算一百分位数

    print(upper_bound)
    #绘制直方图并添加置信区间
    histogram = (
        ggplot(df, aes(x='value', y='..density..')) +
        geom_histogram(bins=30, fill='#0072B2', color='black', alpha=0.8) +
        geom_vline(xintercept=lower_bound, linetype='dashed', color='red') +
        geom_vline(xintercept=upper_bound, linetype='dashed', color='red') +
        geom_density(color='black') +
        labs(title='Xpore', x='Diff mod rate', y='Frequency')+
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
    ggsave(filename=output_path+"oddratio/xpore_0.10.pdf",plot=histogram,dpi=300)
plot_nanocompore()
plot_eligos2()
plot_drummer()
plot_differr()
plot_tombo()
plot_xpore()
