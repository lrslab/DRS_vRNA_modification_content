import pandas as pd
import os
# df=None
# # #
path="SINV/subsets/"
for i in range(1,101):
    file_name=path+'subset'+str(i)+'/DURMMER_result/half2-half1/complete_analysis/NC_001547.1.complete.txt'
    if not os.path.exists(file_name):
        continue
    temp=pd.read_csv(file_name,sep='\t')
    if df is None:
        df=temp
    else:
        df=pd.concat([df, temp], ignore_index=True)
df.to_csv(path+"drummer.bed",sep='\t',index=False)
print(1)
#
df=None
for i in range(1,101):
    file_name=path+'subset'+str(i)+'/eligos2_results/half1_vs_half2_on_gene_baseExt0.txt'
    if not os.path.exists(file_name):
        continue
    temp=pd.read_csv(file_name,sep='\t')
    if df is None:
        df=temp
    else:
        df=pd.concat([df, temp], ignore_index=True)
df.to_csv(path+"eligos2.bed",sep='\t',index=False)
print(1)
df=None
for i in range(1,101):
    file_name=path+'subset'+str(i)+'/xpore_result/diffmod.table'
    temp=pd.read_csv(file_name)
    if df is None:
        df=temp
    else:
        df=pd.concat([df, temp], ignore_index=True)
df.to_csv(path+"xpore.csv",index=False)
print(1)

df=None
for i in range(1,101):

    file_name=path+'subset'+str(i)+'/nanocompore_result/outnanocompore_results.tsv'

    if not os.path.exists(file_name):
        continue
    temp = pd.read_csv(file_name, sep='\t')
    if df is None:
        df=temp
    else:
        df=pd.concat([df, temp], ignore_index=True)
df.to_csv(path+"nanocompore.tsv",index=False,sep='\t')
print(1)