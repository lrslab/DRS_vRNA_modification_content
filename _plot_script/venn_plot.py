from matplotlib_venn import venn2,venn2_circles,venn3
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42
path='subset1_all/c636.txt'
df_all=set(pd.read_csv(path,header=None)[0].values)

path="subset1/c636.txt"
df_772=set(pd.read_csv(path,header=None)[0].values)

path="subset1_500/c636.txt"
df_500=set(pd.read_csv(path,header=None)[0].values)
result_list=[]

my_dpi=300
font2 = { 'size': 8.5}
plt.rc('font', **font2)
plt.figure(figsize=(4,4), dpi=my_dpi)#控制图尺寸的同时，使图高分辨率（高清）显示
out=venn3(subsets = [df_all,df_772,df_500], #绘图数据集
        set_labels = ('All reads', '772x','500x'), #设置组名
        set_colors=("#9790bd","#9ab4db","#a9d7f0"),
        alpha=0.8,#透明度
        normalize_to=1.0,#venn图占据figure的比例，1.0为占满
       )

plt.title("C636")

# temp=df_c636.intersection(df_bhk)
plt.savefig('MIX_subsample/inter_c636.pdf')
plt.show()
print(1)