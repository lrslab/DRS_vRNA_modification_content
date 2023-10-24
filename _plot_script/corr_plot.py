
from plotnine import *
import pandas as pd

from matplotlib import pyplot as plt

plt.rcParams['pdf.fonttype'] = 42

def draw_m6a_net():
       df=pd.read_csv("nanopolish/m6anet_output/data.result.csv")
       df.columns=[0,'Position',1,'Probability',2]
       df=df[['Position','Probability']]
       df.columns=['Position','C636']
       df_ivt=pd.read_csv("m6anet_output/data.result.csv")
       df_ivt.columns=[0,'Position',1,'Probability',2]
       df_ivt=df_ivt[['Position','Probability']]
       df_ivt.columns=['Position','IVT']

       final_df=df.merge(df_ivt,on='Position',how='inner')

       df = pd.DataFrame({'x': [0, 1], 'y': [0, 1]})

       plot = ggplot(final_df, aes(x='C636', y='IVT')) +\
              geom_point(alpha = 0.5,color='none',fill='black',size=1) \
              + theme_bw() \
              + geom_line(df,aes(x='x', y='y'),color='red',alpha = 0.9,linetype='dashed')\
              + ylim(0, 1) \
              + xlim(0, 1)\
              + theme(
                      figure_size=(4,4),
                      axis_text=element_text(size=13,
                                             family="Arial"),
                      axis_title=element_text(size=13,
                                              family="Arial"),
                      panel_grid_minor=element_blank(),
                      title=element_text(size=13,
                                         family="Arial")
                      )\
              +labs(title='M6anet')

       ggsave(filename="machine_learning_plot/C636_m6anet.pdf",plot=plot,dpi=300)
       print(plot)


def draw_nanom6A():
       df=pd.read_csv("ratio.bed",sep='\t',header=None)
       df.columns=[0,'Position',1,2,3,4,5,'Probability']
       df = df[['Position', 'Probability']]
       df.columns = ['Position', 'C636']
       df_ivt=pd.read_csv("m6A_ratio.bed",sep='\t',header=None)
       df_ivt.columns=[0,'Position',1,2,3,4,5,'Probability']
       df_ivt = df_ivt[['Position', 'Probability']]
       df_ivt.columns = ['Position', 'IVT']

       final_df = df.merge(df_ivt, on='Position', how='inner')

       df = pd.DataFrame({'x': [0, 1], 'y': [0, 1]})

       plot = ggplot(final_df, aes(x='C636', y='IVT')) + \
              geom_point(alpha=0.5, color='none', fill='black',size=1) \
              + theme_bw() \
              + geom_line(df, aes(x='x', y='y'), color='red', alpha=0.9, linetype='dashed') \
              + ylim(0, 1) \
              + xlim(0, 1) \
              + theme(
              figure_size=(4, 4),
              axis_text=element_text(size=13,
                                     family="Arial"),
              axis_title=element_text(size=13,
                                      family="Arial"),
              panel_grid_minor=element_blank(),
              title=element_text(size=13,
                                 family="Arial")
       ) \
              + labs(title='Nanom6A')

       ggsave(filename="machine_learning_plot/C636_nanom6A.pdf", plot=plot,
              dpi=300)
       print(plot)

def draw_tombo():
       def extract_tombo_denovo_result(path):
              file = open(path, 'r')
              result_list = []
              for line in file:
                     if line[0] == '>':
                            temp_list = line[:-1].split(':')
                            item = [int(temp_list[1]) - 1, float(temp_list[-1])]
                     else:
                            if line[0] == 'A':
                                   result_list.append(item)
              return pd.DataFrame(result_list)

       denovo_path = "tombo_results.significant_regions.fasta"
       df = extract_tombo_denovo_result(denovo_path)
       df.columns = ['Position', 'C636']
       denovo_path = "tombo_results.significant_regions.fasta"
       df_ivt = extract_tombo_denovo_result(denovo_path)
       df_ivt.columns = ['Position', 'IVT']
       final_df = df.merge(df_ivt, on='Position', how='inner')

       df = pd.DataFrame({'x': [0, 1], 'y': [0, 1]})

       plot = ggplot(final_df, aes(x='C636', y='IVT')) + \
              geom_point(alpha=0.5, color='none', fill='black',size=1) \
              + theme_bw() \
              + geom_line(df, aes(x='x', y='y'), color='red', alpha=0.9, linetype='dashed') \
              + ylim(0, 1) \
              + xlim(0, 1) \
              + theme(
              figure_size=(4, 4),
              axis_text=element_text(size=13,
                                     family="Arial"),
              axis_title=element_text(size=13,
                                      family="Arial"),
              panel_grid_minor=element_blank(),
              title=element_text(size=13,
                                 family="Arial")
       ) \
              + labs(title='Tombo')

       ggsave(filename="machine_learning_plot/C636_tombo.pdf", plot=plot,
              dpi=300)
       print(plot)
       print(1)
draw_m6a_net()
draw_nanom6A()
draw_tombo()
# correlation