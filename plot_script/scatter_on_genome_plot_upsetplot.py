
import numpy as np
from plotnine import *
import pandas as pd
from collections import  OrderedDict

from matplotlib import pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
fasta_path="SINV_Toto1101.fa"
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
    df=df[df["oddR"]>2.34]
    df = df[df["-logFDR"] > 1.3]
    df[3] = df.apply(lambda x: fasta[x["start"]], axis=1)
    df = df[df[3] == 'C']
    return df

ref_dict={
    "A":"T","T":"A","C":"G","G":"C"
}
fasta=list(read_fasta_to_dic(fasta_path).values())[0]
def extract_eligos2(eligos2_path):
    df=pd.read_csv(eligos2_path,sep='\t')
    df['LogOddRatio'] = np.log2(df["oddR"].values)
    df['-logPvalue'] = np.log(df['pval']) * (-1)
    df=df[(df['-logPvalue'] >3) & (df['LogOddRatio'] >0.83)]
    df=df[df["ref"]=='C']
    df['start']=df['start_loc']
    # df = df.replace(np.inf, np.nan).dropna(subset=['LogOddRatio'])
    return df
def test_genome(fasta,df):
    ref_dict={"A":"T","T":"A","C":"G","G":"C"}
    df["ref"]=df.apply(lambda x: True if fasta[x["position"]]==x['ref']  else False,axis=1)
    return df

def extract_xpore(path):
    df_original = pd.read_csv(path)
    genome_fasta=read_fasta_to_dic(fasta_path).values()
    genome_fasta=list(genome_fasta)[0]
    df_original["-logPvalue"]=df_original["pval_IVT_vs_WT"].apply(lambda x:-np.log10(x))

    df_original["end_position"]=df_original["position"]+1

    df_original["ref"] = df_original.apply(lambda x: x["kmer"][2], axis=1)
    df_original=df_original[df_original['ref']=='C']

    df_original=test_genome(genome_fasta,df_original)
    df_original=df_original[df_original["ref"]==True]
    df_original.drop("ref",inplace=True,axis=1)
    df_original['diff_mod_rate'] = df_original["diff_mod_rate_IVT_vs_WT"]
    df_original = df_original[df_original['-logPvalue'] > 4]
    df_original = df_original[(df_original['diff_mod_rate'] <-0.08)|(df_original['diff_mod_rate'] >0.08)]
    df_original["start"]=df_original["position"]
    return df_original

def read_nanocompore(path):
    df_original = pd.read_csv(path,sep='\t')
    genome_fasta=read_fasta_to_dic(fasta_path).values()
    genome_fasta=list(genome_fasta)[0]
    # df_original["GMM_logit_pvalue"]=df_original["GMM_logit_pvalue"].apply(lambda x:-np.log10(x))
    df_original["ref"]=df_original.apply(lambda x:x["ref_kmer"][0],axis=1)
    # df_result["."]="."
    # df_result=df_result[["id","position","end_position",'.',"pval_IVT_vs_L",'strand']]
    df_original['-logPvalue'] = np.log10(df_original['GMM_logit_pvalue']) * (-1)
    df_original=df_original[df_original['-logPvalue']>2]
    df_original = df_original.replace("NC", np.nan).dropna(subset=['Logit_LOR'])
    df_original['LogOddRatio'] = df_original['Logit_LOR'].astype(float)
    df_original = df_original[(df_original['LogOddRatio'] <-0.5)|(df_original['LogOddRatio'] >0.5)]
    df_original=df_original[df_original['ref']=='C']
    # df_original=test_genome(genome_fasta,df_original)
    # df_original=df_original[df_original["ref"]==True]
    df_original.drop("ref",inplace=True,axis=1)
    df_original['start']=df_original['pos']
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
    all_df.columns=['C', 0, 'A1', 'A2', '-log10(P_value)', 3 ,2, 'difference']
    all_df=all_df[['C', 0, 'A1', 'A2', '-log10(P_value)', 3 , 'difference',2]]

    all_df['sig'] = 'normal'
    # sns.displot(data=all_df, x="-log10(P_value)", kind="kde")
    # # plt.xscale('log')
    # plt.show()
    all_df = all_df[all_df[3] == 'C']
    all_df.loc[(all_df.difference>  0.05 )&(all_df["-log10(P_value)"] >2),'sig'] = 'up'
    all_df.loc[(all_df.difference< -0.05)&(all_df["-log10(P_value)"] >2),'sig'] = 'down'
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

    all_df=all_df[(all_df["sig"]=="up")|(all_df["sig"]=="down")]
    # all_df = all_df[all_df["-log10(P_value)"] >2]
    all_df.columns = ["chro", 'start', 'end', '.', '-logPvalue', '..', 'difference', 'strand',
                          'sig']
    return all_df

def extract_drummer(drummer_path):
    df= pd.read_csv(drummer_path, sep='\t')
    df['-logPvalue'] = np.log10(df['p_values_OR_adj']) * (-1)
    df=df[df['p_values_OR_adj']<0.05]
    df['LogOddRatio'] = np.log2(df["odds_ratio"].values)
    df=df[df["LogOddRatio"]>0.93]
    df=df[df['ref_ctrl']=='C']
    df['start']=df['pos_ctrl']-1
    return df

original_path="result/"
differr_path=original_path+"differr.bed"
differr_path=extract_differr(differr_path)

eligos2_path=original_path+"eligos2_results/bhk_vs_ivt_on_gene_baseExt0.txt"
eligos2_path=extract_eligos2(eligos2_path)
#
nanocompore_path= original_path+"nanocompore_result/outnanocompore_results.tsv"
nanocompore_path=read_nanocompore(nanocompore_path)
#
tombo_path=original_path+"sample"
tombo_path=extract_tombo(tombo_path)
#
xpore_path=original_path+"xpore_result/diffmod.table"
xpore_path=extract_xpore(xpore_path)

drummer_path=original_path+"DURMMER_result/bhk-ivt/complete_analysis/NC_001547.1.complete.txt"
drummer_path=extract_drummer(drummer_path)
output_path="plot/"
def plot_differr(differr_path):
    cutoff = 2
    p = ggplot(aes(x='start', y='LogOddRatio',color="-logFDR"), differr_path) \
        + geom_point()\
        + theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (6, 3),
                legend_text=element_text(family="Arial"),
                legend_title=element_text(family="Arial"),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13, family="Arial"),
                axis_title=element_text(size=13, family="Arial"),
                title=element_text(size=13,
                                   family="Arial")
                )\
        + coord_cartesian(xlim =(0, 11703))\
        + scale_x_continuous(breaks = (0,2000,4000,6000,8000,10000,11703)) \
        + geom_hline(yintercept=2.34, linetype='dashed', color='black')\
        + scale_color_gradient(low="#F2F2F2", high="red")\
        +labs(title='Differr', x='Position', y='Log2(Odds Ratio)')
    ggsave(filename=output_path+"differr.pdf",plot=p,dpi=300)
    print(p)
plot_differr(differr_path)

def plot_eligos2(eligos2_path):
    cutoff = 2
    p = ggplot(aes(x='start', y='LogOddRatio',color="-logPvalue"), eligos2_path) \
        + geom_point()\
        + theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (6, 3),
                legend_text=element_text(family="Arial"),
                legend_title=element_text(family="Arial"),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13, family="Arial"),
                axis_title=element_text(size=13, family="Arial"),
                title=element_text(size=13,
                                   family="Arial")
                )\
        + coord_cartesian(xlim =(0, 11703))\
        + scale_x_continuous(breaks = (0,2000,4000,6000,8000,10000,11703)) \
        + geom_hline(yintercept=0.83, linetype='dashed', color='black')\
        + scale_color_gradient(low="#F2F2F2", high="red")\
        +labs(title='ELIGOS2', x='Position', y='Log2(Odds Ratio)')
    ggsave(filename=output_path+"eligos2.pdf",plot=p,dpi=300)
    print(p)
print(2)
plot_eligos2(eligos2_path)

def plot_xpore(xpore_path):
    cutoff = 2
    p = ggplot(aes(x='start', y='diff_mod_rate',color="-logPvalue"), xpore_path) \
        + geom_point()\
        + theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (6, 3),
                legend_text=element_text(family="Arial"),
                legend_title=element_text(family="Arial"),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13, family="Arial"),
                axis_title=element_text(size=13, family="Arial"),
                title=element_text(size=13,
                                   family="Arial")
                )\
        + coord_cartesian(xlim =(0, 11703))\
        + scale_x_continuous(breaks = (0,2000,4000,6000,8000,10000,11703)) \
        + scale_color_gradient(low="#F2F2F2", high="red") \
        + geom_hline(yintercept=-0.08, linetype='dashed', color='black') \
        + geom_hline(yintercept=0.08, linetype='dashed', color='black')\
        +labs(title='Xpore', x='Position', y='Differential modification rate')
    ggsave(filename=output_path+"xpore.pdf",plot=p,dpi=300)
    print(p)
plot_xpore(xpore_path)
print(2)
def plot_tombo(tombo_path):
    cutoff = 2
    p = ggplot(aes(x='start', y='difference',color="-logPvalue"), tombo_path) \
        + geom_point()\
        + theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (6, 3),
                legend_text=element_text(family="Arial"),
                legend_title=element_text(family="Arial"),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13, family="Arial"),
                axis_title=element_text(size=13, family="Arial"),
                title=element_text(size=13,
                                   family="Arial")
                )\
        + coord_cartesian(xlim =(0, 11703))\
        + scale_x_continuous(breaks = (0,2000,4000,6000,8000,10000,11703)) \
        + geom_hline(yintercept=-0.05, linetype='dashed', color='black') \
        + geom_hline(yintercept=0.05, linetype='dashed', color='black') \
        + scale_color_gradient(low="#F2F2F2", high="red")\
        +labs(title='Tombo', x='Position', y='Difference')
    ggsave(filename=output_path+"tombo.pdf",plot=p,dpi=300)
    print(p)
plot_tombo(tombo_path)
def plot_nanocompore(nanocompore_path):
    cutoff = 2
    p = ggplot(aes(x='start', y='LogOddRatio',color="-logPvalue"), nanocompore_path) \
        + geom_point()\
        + theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                figure_size = (6, 3),
                legend_text=element_text(family="Arial"),
                legend_title=element_text(family="Arial"),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13, family="Arial"),
                axis_title=element_text(size=13, family="Arial"),
                title=element_text(size=13,
                                   family="Arial")
                )\
        + coord_cartesian(xlim =(0, 11703)) \
        + geom_hline(yintercept=0.5, linetype='dashed', color='black') \
        + geom_hline(yintercept=-0.5, linetype='dashed', color='black')\
        + scale_x_continuous(breaks = (0,2000,4000,6000,8000,10000,11703)) \
        + scale_color_gradient(low="#F2F2F2", high="red")\
        +labs(title='Nanocompore', x='Position', y='Log2(Odds Ratio)')
    # ggsave(filename=output_path+"nanocompore.svg",plot=p,dpi=300)
    ggsave(filename=output_path+"nanocompore.pdf", plot=p, dpi=300)
    # ggplot.save(output_path+"nanocompore.ai")
    print(p)
plot_nanocompore(nanocompore_path)

def plot_drummer(drummer_path):
    cutoff = 2
    p = ggplot(aes(x='start', y='LogOddRatio',color="-logPvalue"), drummer_path) \
        + geom_point()\
        + theme_bw()\
        + theme(legend_key_height=10,
                legend_key_width=5,
                legend_text=element_text(family="Arial"),
                legend_title=element_text(family="Arial"),
                figure_size = (6, 3),
                panel_grid_minor=element_blank(),
                axis_text=element_text(size=13, family="Arial"),
                axis_title=element_text(size=13, family="Arial"),
                title=element_text(size=13,
                                   family="Arial")

                )\
        + coord_cartesian(xlim =(0, 11703))\
        + scale_x_continuous(breaks = (0,2000,4000,6000,8000,10000,11703)) \
        + geom_hline(yintercept=0.93, linetype='dashed', color='black')\
        + scale_color_gradient(low="#F2F2F2", high="red")\
        +labs(title='DRUMMER', x='Position', y='Log2(Odds Ratio)')
    ggsave(filename=output_path+"drummer.pdf",plot=p,dpi=300)
    print(p)
plot_drummer(drummer_path)
from upsetplot import from_contents
from upsetplot import UpSet
from matplotlib import pyplot as plt
plt.rcParams['font.family'] = 'Arial'
result_dict=OrderedDict()
result_dict["Differr"]=differr_path["start"].values
result_dict["DRUMMER"]=drummer_path["start"].values
result_dict["ELIGOS2"]=eligos2_path["start_loc"].values
result_dict["Nanocompore"]=nanocompore_path["pos"].values
result_dict["Tombo"]=tombo_path["start"].values
result_dict["Xpore"]=xpore_path["position"].values

animals = from_contents(result_dict)
UpSet(animals, subset_size='count',show_counts=True,sort_categories_by='-input').plot()


plt.savefig(output_path+'upsetplot_C.pdf',dpi=300)
plt.show()
tombo=set(tombo_path["start"].values)
xpore=set(xpore_path["position"].values)
temp=tombo.intersection(xpore)
print(1)