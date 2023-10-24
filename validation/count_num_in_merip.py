import pandas as pd

site_df = pd.read_csv('data/prediction',sep='\t',header=None)
inter_df = pd.read_csv('data/merip.bed',sep='\t',header=None)
# inter_df=inter_df[inter_df[6]=='Vero']
i=1
for index,line in site_df.iterrows():
    pos = int(line[1])
    temp = inter_df[(inter_df[1]<=pos)&(inter_df[2]>=pos)]
    if temp.shape[0]>=1:
        i=i+1
print(i)