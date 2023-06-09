
import numpy as np
from plotnine import *
import pandas as pd
from collections import  OrderedDict
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
fasta_path="/t1/zhguo/Data/Virus_data/SINV/SINV_Toto1101.fa"
def read_fasta_to_dic(filename):
    """
    function used to parser small fasta
    still effective for genome level file
    """
    fa_dic = OrderedDict()

    with open(filename, "r") as f:
        for n, line in enumerate(f.readlines()):
            if line.startswith(">"):
                if n > 0:
                    fa_dic[short_name] = "".join(seq_l)  # store previous one

                full_name = line.strip().replace(">", "")
                short_name = full_name.split(" ")[0]
                seq_l = []
            else:  # collect the seq lines
                if len(line) > 8:  # min for fasta file is usually larger than 8
                    seq_line1 = line.strip()
                    seq_l.append(seq_line1)

        fa_dic[short_name] = "".join(seq_l)  # store the last one
    return fa_dic

def extract_differr(path):
    df=pd.read_csv(path,sep='\t',header=None)
    df = df.iloc[:, 0:10]
    df.columns = ["chro", 'start', 'end', 'ref', '.', 'strand', 'oddR', 1, 2, "-logFDR"]
    df['LogOddRatio'] = df["oddR"]
    # df=df[df["oddR"]>1]
    # df = df[df["-logFDR"] > 1.3]
    df[3] = df.apply(lambda x: fasta[x["start"]], axis=1)
    # df = df[df[3] == 'A']
    df["test"]=(df["-logFDR"] > 1.3)& (df["oddR"]>1)
    return df

ref_dict={
    "A":"T","T":"A","C":"G","G":"C"
}
fasta=list(read_fasta_to_dic(fasta_path).values())[0]
def extract_eligos2(eligos2_path):
    df=pd.read_csv(eligos2_path,sep='\t')
    df['LogOddRatio'] = np.log2(df["oddR"].values)
    df['-logPval'] = np.log(df['pval']) * (-1)
    # df=df[(df['-logPval'] >=3) & (df['LogOddRatio'] >= 1.2)]
    # df=df[df["ref"]=='A']
    df['start']=df['start_loc']
    # df = df.replace(np.inf, np.nan).dropna(subset=['LogOddRatio'])
    df["test"] = (df["-logPval"] > 3) & (df["oddR"] > 1.2)
    return df
def test_genome(fasta,df):
    ref_dict={"A":"T","T":"A","C":"G","G":"C"}
    df["ref"]=df.apply(lambda x: True if fasta[x["position"]]==x['ref']  else False,axis=1)
    return df

def extract_xpore(path):
    df_original = pd.read_csv(path)
    genome_fasta=read_fasta_to_dic(fasta_path).values()
    genome_fasta=list(genome_fasta)[0]
    df_original["-logPval"]=df_original["pval_IVT_vs_WT"].apply(lambda x:-np.log10(x))

    df_original["end_position"]=df_original["position"]+1

    df_original["ref"] = df_original.apply(lambda x: x["kmer"][2], axis=1)
    # df_original=df_original[df_original['ref']=='A']

    df_original=test_genome(genome_fasta,df_original)
    df_original=df_original[df_original["ref"]==True]
    df_original.drop("ref",inplace=True,axis=1)
    df_original['diff_mod_rate'] = df_original["diff_mod_rate_IVT_vs_WT"]
    # df_original = df_original[df_original['-logPval'] > 2]
    # df_original = df_original[(df_original['diff_mod_rate'] <-0.08)|(df_original['diff_mod_rate'] >0.05)]
    df_original["start"]=df_original["position"]
    df_original["test"] = df_original["-logPval"] > 3
    return df_original

def read_nanocompore(path):
    df_original = pd.read_csv(path,sep='\t')
    genome_fasta=read_fasta_to_dic(fasta_path).values()
    genome_fasta=list(genome_fasta)[0]
    # df_original["GMM_logit_pvalue"]=df_original["GMM_logit_pvalue"].apply(lambda x:-np.log10(x))
    df_original["ref"]=df_original.apply(lambda x:x["ref_kmer"][0],axis=1)
    # df_result["."]="."
    # df_result=df_result[["id","position","end_position",'.',"pval_IVT_vs_L",'strand']]
    df_original['-logPval'] = np.log10(df_original['GMM_logit_pvalue']) * (-1)
    # df_original=df_original[df_original['-logPval']>2]
    df_original = df_original.replace("NC", np.nan).dropna(subset=['Logit_LOR'])
    df_original['LogOddRatio'] = df_original['Logit_LOR'].astype(float)
    # df_original = df_original[(df_original['LogOddRatio'] <-0.406)|(df_original['LogOddRatio'] >0.406)]
    # df_original=df_original[df_original['ref']=='A']
    # df_original=test_genome(genome_fasta,df_original)
    # df_original=df_original[df_original["ref"]==True]
    df_original.drop("ref",inplace=True,axis=1)
    df_original['start']=df_original['pos']
    df_original["test"] = (df_original["-logPval"] > 3) & ((df_original["LogOddRatio"] > 0.5)|(df_original["LogOddRatio"] <- 0.5))
    return  df_original


def extract_tombo(ab_path):
    tombo_result_minus_path= ab_path + ".level_samp_comp_detect.statistic.minus.wig"
    tombo_result_plus_path= ab_path +".level_samp_comp_detect.statistic.plus.wig"
    tombo_difference_minus_path= ab_path + ".level_samp_comp_detect.difference.minus.wig"
    tombo_difference_plus_path= ab_path + ".level_samp_comp_detect.difference.plus.wig"
    fasta=list(read_fasta_to_dic(fasta_path).values())[0]
    ref_dict={
        "A":"T",
        "T":"A",
        "C":"G",
        "G":"C"
    }
    def read_tombo_results(path,strand):
        df = pd.read_csv(path,sep=' ',skiprows=2, header=None)
        # df = df[df[1]>=0.2]
        df[2]=df[0].apply(lambda x: fasta[x] if strand=='+' else ref_dict[fasta[x]])
        df.insert(loc=0, column="A", value=np.repeat([["NC_000913.3"]],df.shape[0]))
        df.insert(loc=2, column="A1", value=df[0].values+1)
        df.insert(loc=3, column="A2", value=np.repeat([["."]], df.shape[0]))
        # sns.displot(df, x=1)
        # plt.show()
        if path.find("plus") != -1:
            # df = df[df[2] == "A"]
            df[3]='+'
        else:
            # df = df[df[2] == "T"]
            df[3]='-'
        return df

    plus_df = read_tombo_results(tombo_difference_plus_path,'+')
    try:
        minus_df = read_tombo_results(tombo_difference_minus_path,'-')
    except Exception:
        minus_df=None
    df_new = pd.concat([plus_df,minus_df], ignore_index=True)


    plus_df = read_tombo_results(tombo_result_plus_path,'+')
    try:
        minus_df = read_tombo_results(tombo_result_minus_path,'-')
    except Exception:
        minus_df=None
    df = pd.concat([plus_df,minus_df], ignore_index=True)


    all_df=pd.merge(df,df_new,how='inner',on=["A",0,"A1","A2",3,2])
    all_df.columns=['A', 0, 'A1', 'A2', '-log10(P_value)', 3 ,2, 'difference']
    all_df=all_df[['A', 0, 'A1', 'A2', '-log10(P_value)', 3 , 'difference',2]]

    all_df['sig'] = 'normal'
    # sns.displot(data=all_df, x="-log10(P_value)", kind="kde")
    # # plt.xscale('log')
    # plt.show()
    # all_df = all_df[all_df[3] == 'A']
    # all_df.loc[(all_df.difference>  0.047 )&(all_df["-log10(P_value)"] >2),'sig'] = 'up'
    # all_df.loc[(all_df.difference< -0.047)&(all_df["-log10(P_value)"] >2),'sig'] = 'down'
    # base_plot = ggplot(all_df) +\
    #             geom_point(aes(x='difference', y='-log10(P_value)',color="sig"),alpha=0.3) + \
    #             scale_color_manual(values ={"up":"red","down":"blue","normal":"black"}) +\
    #             theme_bw() + \
    #             theme(figure_size=(6, 6),
    #                   legend_box_spacing=0.1,
    #                   legend_background=element_blank(),
    #                   legend_direction='vert',
    #                   legend_title=element_blank()) +\
    #             coord_cartesian(xlim=(-1, 1), expand=True)
    # print(base_plot)

    # all_df=all_df[(all_df["sig"]=="up")|(all_df["sig"]=="down")]
    # all_df = all_df[all_df["-log10(P_value)"] >2]

    all_df.columns = ["chro", 'start', 'end', '.', '-logPval', '..', 'difference', 'strand',
                          'sig']
    all_df['test']=all_df['-logPval']>2
    return all_df

def extract_drummer(drummer_path):
    df= pd.read_csv(drummer_path, sep='\t')
    df['-logPvalue'] = np.log10(df['p_values_OR_adj']) * (-1)
    df["test"] = (df['p_values_OR_adj']<0.05) & (df["odds_ratio"]>1.5)
    df['LogOddRatio'] = np.log2(df["odds_ratio"].values)
    return df

original_path="/t1/zhguo/Data/Virus_data/SINV/IVT/different_coverage/20_reads/"
original_path='/t1/zhguo/Data/Virus_data/SINV/IVT/100_subsample_gRNA/subset1/'
differr_path=original_path+"differr.bed"
differr_path=extract_differr(differr_path)

eligos2_path=original_path+"eligos2_results/ivt1_vs_ivt2_on_gene_baseExt0.txt"
eligos2_path=extract_eligos2(eligos2_path)
#
nanocompore_path= original_path+"nanocompore_results/outnanocompore_results.tsv"
nanocompore_path=read_nanocompore(nanocompore_path)
#
tombo_path=original_path+"sample"
tombo_path=extract_tombo(tombo_path)
#
xpore_path=original_path+"xpore_result/diffmod.table"
xpore_path=extract_xpore(xpore_path)

drummer_path=original_path+"DURMMER_result/ivt2-ivt1/complete_analysis/NC_001547.1.complete.txt"
drummer_path=extract_drummer(drummer_path)

def plot_drummer(drummer_path):
    cutoff = 2
    p = ggplot(aes( x='LogOddRatio',y="-logPvalue",fill='test'), drummer_path) \
        + geom_point(alpha=0.5,color='none')\
        + theme_bw()\
        + theme(
                figure_size = (4,4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                legend_position='none',
                title = element_text(size=13,
                                     family="Arial")
                ) \
        + xlim(-4, 4) \
        + ylim(0, 2) \
        + geom_vline(xintercept=np.log2(1.5), linetype='dashed', color='black') \
        + geom_hline(yintercept=-np.log10(0.05), linetype='dashed', color='black') \
        + scale_fill_manual(values=("#A9A9A9","red"))\
        +labs(title='DRUMMER' ,y="-logPvalue", x='LogOddRatio')
    ggsave(filename="/t1/zhguo/Data/Virus_data/SINV/IVT/plot_result/point_cutoff/drummer.pdf",plot=p,dpi=300)
    print(p)
plot_drummer(drummer_path)

def plot_differr(differr_path):
    cutoff = 2
    p = ggplot(aes(x='LogOddRatio', y='-logFDR',fill='test'), differr_path) \
        + geom_point(alpha=0.5,color='none')\
        + theme_bw()\
        + theme(
                figure_size = (4, 4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                title=element_text(size=13,
                                   family="Arial"),
                legend_position='none'
                ) \
        + xlim(-16, 16) \
        + geom_vline(xintercept=1, linetype='dashed', color='black') \
        + geom_hline(yintercept=1.3, linetype='dashed', color='black') \
        + scale_fill_manual(values=("#A9A9A9","red"))\
        +labs(title='Differr', y='-LogFDR', x='LogOddRatio')
    ggsave(filename="/t1/zhguo/Data/Virus_data/SINV/IVT/different_coverage/20_reads/plot/differr.pdf",plot=p,dpi=300)
    print(p)
plot_differr(differr_path)
def plot_eligos2(eligos2_path):
    cutoff = 2
    p = ggplot(aes(y='-logPval', x='LogOddRatio',fill="test"), eligos2_path) \
        + geom_point(alpha=0.5,color='none')\
        + theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (4, 4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                title=element_text(size=13,
                                   family="Arial"),
                legend_position='none'
                )\
        + xlim(-4,4)\
        + geom_vline(xintercept=np.log2(1.2), linetype='dashed', color='black') \
        + geom_hline(yintercept=3, linetype='dashed', color='black') \
        + scale_fill_manual(values=("#A9A9A9", "red")) \
        +labs(title='ELIGOS2', y='-LogPvalue', x='LogOddRatio')
    ggsave(filename="/t1/zhguo/Data/Virus_data/SINV/IVT/different_coverage/20_reads/plot/eligos2.pdf",plot=p,dpi=300)
    print(p)
print(2)
plot_eligos2(eligos2_path)

def plot_xpore(xpore_path):
    cutoff = 2
    p = ggplot(aes(y='-logPval', x='diff_mod_rate',fill="test"), xpore_path) \
        + geom_point(alpha=0.5,color='none')\
        + theme_bw()\
        + theme(
                figure_size = (4, 4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                title=element_text(size=13,
                                   family="Arial"),
                legend_position='none'
                ) \
        + xlim(-0.2, 0.2) \
        + scale_fill_manual(values=("#A9A9A9","red")) \
        + geom_hline(yintercept=3, linetype='dashed', color='black') \
        +labs(title='Xpore', y='-LogPvalue', x='Differential modification rate')
    ggsave(filename="/t1/zhguo/Data/Virus_data/SINV/IVT/different_coverage/20_reads/plot/xpore.pdf",plot=p,dpi=300)
    print(p)
plot_xpore(xpore_path)
print(2)
def plot_tombo(tombo_path):
    cutoff = 2
    p = ggplot(aes(y='-logPval', x='difference',fill="test"), tombo_path) \
        + geom_point(alpha=0.5,color='none')\
        + theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (4, 4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                title=element_text(size=13,
                                   family="Arial"),
                legend_position='none'
                ) \
        + xlim(-0.2, 0.2) \
        + scale_fill_manual(values=("#A9A9A9", "red")) \
        + geom_hline(yintercept=2, linetype='dashed', color='black') \
        +labs(title='Tombo', x='Difference', y='-LogPvalue')
    ggsave(filename="/t1/zhguo/Data/Virus_data/SINV/IVT/different_coverage/20_reads/plot/tombo.pdf",plot=p,dpi=300)
    print(p)
plot_tombo(tombo_path)
def plot_nanocompore(nanocompore_path):
    cutoff = 2
    p = ggplot(aes(x='LogOddRatio',y="-logPval",fill='test'), nanocompore_path) \
        + geom_point(alpha=0.5,color='none')\
        + theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (4, 4),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13,
                                              family="Arial"),
                axis_title=element_text(size=13,
                                              family="Arial"),
                title=element_text(size=13,
                                   family="Arial"),
                legend_position='none'
                )\
        + ylim(0,5)\
        + geom_vline(xintercept=0.5, linetype='dashed', color='black') \
        + geom_vline(xintercept=-0.5, linetype='dashed', color='black') \
        + geom_hline(yintercept=2, linetype='dashed', color='black') \
        +labs(title='Nanocompore', y='-LogPvalue', x='LogOddRatio') \
        + scale_fill_manual(values=("#A9A9A9"))
    ggsave(filename="/t1/zhguo/Data/Virus_data/SINV/IVT/different_coverage/20_reads/plot/nanocompore.pdf",plot=p,dpi=300)
    print(p)
plot_nanocompore(nanocompore_path)
from upsetplot import from_contents
from upsetplot import UpSet
result_dict=OrderedDict()
result_dict["differr"]=differr_path["start"].values
result_dict["eligos2"]=eligos2_path["start_loc"].values
result_dict["nanocompore"]=nanocompore_path["pos"].values
result_dict["tombo"]=tombo_path["start"].values
result_dict["xpore"]=xpore_path["position"].values

animals = from_contents(result_dict)
UpSet(animals, subset_size='count',show_counts=True).plot()

from matplotlib import pyplot as plt
plt.show()
print(1)