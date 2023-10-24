import pandas as pd
from plotnine import *
import pandas as pd

from matplotlib import pyplot as plt

plt.rcParams['pdf.fonttype'] = 42

df = pd.read_csv('data/merip.bed',sep='\t',header=None)

df_site = pd.read_csv('data/prediction',sep='\t',header=None)
df_site[6] = 'ONT'
df = pd.concat([df,df_site],axis=0)
df.columns=[0,'x_start','x_end',3,4,5,'group']
df['y_start'] = df['group'].apply(lambda x: 3.25 if x=='Caco-2' else 2.25 if x=='Vero' else 1.25)
df['y_end'] = df['group'].apply(lambda x: 2.75 if x=='Caco-2' else 1.75 if x=='Vero' else 0.75)

p=ggplot(df,aes(xmin='x_start', xmax='x_end', ymin='y_start', ymax='y_end',fill='group')) + \
    geom_rect() + \
    theme_bw() +\
    xlim(0,29900)+\
    labs(title='Multiple Rectangles', x='X', y='Y')+\
    theme(
      figure_size=(20,3),
      axis_text=element_text(size=13,
                             family="Arial"),
      axis_title=element_text(size=13,
                              family="Arial"),
      panel_grid_minor=element_blank(),
      title=element_text(size=13,
                         family="Arial"),
      legend_position='none'
      )\

print(p)
p.save(filename="sars.pdf",dpi=300)
print(1)