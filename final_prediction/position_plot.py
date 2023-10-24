import pandas as pd
import numpy as np
import plotnine as p9
from matplotlib import pyplot as plt

plt.rcParams['pdf.fonttype'] = 42

c636 = pd.read_csv('data/sgRNA/c636_intersection.bed', sep='\t')[['start','ref']]
c636['group1'] = 'C6/36'
bhk = pd.read_csv('data/sgRNA/bhk_intersection.bed', sep='\t')[['start','ref']]
bhk['group2'] = 'BHK-21'
final = pd.merge(c636, bhk, on=['start','ref'], how='outer')
final['group'] = final.apply(
    lambda x: 'inter' if pd.notna(x['group1']) and pd.notna(x['group2']) else x['group1'] if pd.notna(x['group1']) else
    x['group2'], axis=1)
final['y'] = final.apply(
    lambda x: 1 if x['group'] == 'BHK-21' or x['group'] == 'inter' else 0, axis=1)
final['yend'] = final.apply(
    lambda x: -1 if x['group'] == 'C6/36' or x['group'] == 'inter' else 0, axis=1)

p9plot = p9.ggplot(final, p9.aes(x='start', y='y', xend='start', yend='yend',color='ref')) + \
    p9.geom_segment() + \
    p9.geom_vline(xintercept=11703)+\
    p9.theme_bw()+ \
    p9.xlim(0, 11703) + \
    p9.scale_y_continuous(limits=(-3, 3)) + \
    p9.scale_colour_manual(values={
                          "A":'#00CC00',
                          "T":'#CC0000',
                          "G":'#FFB300',
                          "C":"#0000CC"})+\
    p9.theme(
        figure_size=(10, 3),
        panel_grid_minor=p9.element_blank(),
        axis_text=p9.element_text(size=12,family='Arial'),
        axis_title=p9.element_text(size=12,family='Arial'),
        title=p9.element_text(size=12,family='Arial'),

        legend_position='none',


)
print(p9plot)
p9plot.save('sgRNA.pdf',dpi=300)