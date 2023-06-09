import pandas as pd
import os
import numpy as np
from tqdm import tqdm

def select_sub_reads(line):
    result = line[0].split("/")[-1].split(".")[0]
    return result

def copy2new_path(line):
    global pbar
    fold_name = line[0].split('/')[-2]
    if not os.path.exists(fold_name):
        os.mkdir(fold_name)
    os.system('cp ' + line[0] + ' ' + fold_name)
    pbar.update(1)

# original fast5(single format)
fast5_path = "SINV/IVT/single"
# output fast5 path
results_path = 'SINV/IVT/subsets/'
# read_list
df = pd.read_csv("SINV_IVT_FL.csv")

def transfer_fast5(total_fl, df ,path):
    global pbar
    total_fl.columns = [0, 'readname']
    df = df[['readname']]
    df[2] = '.'
    method='inner'
    if method =='inner':
        temp = pd.merge(total_fl, df, on="readname", how='inner')
    else:
        temp = pd.merge(total_fl, df, on="readname", how='outer')
        temp=temp[temp[2].isna()]
    os.chdir(path)
    cmd='mkdir single'
    os.system(cmd)
    os.chdir(path+'/single')
    print(os.getcwd())
    os.system("rm -rf *")
    pbar = tqdm(total=temp.shape[0], position=0, leave=True)
    temp.apply(copy2new_path, axis=1)

#READ FAST5 FILE LIST

os.system("find " + fast5_path + " -name \"*.fast5\" >" + "files.txt")
fast5_file = 'files.txt'
print("Generated fast5 list file: ", fast5_file)
total_fl = []
for i in open(fast5_file, "r"):
    total_fl.append(i.rstrip())
total_fl = pd.DataFrame(np.array(total_fl))

total_fl[1] = total_fl.apply(select_sub_reads, axis=1)

for i in range(1,100):
    os.chdir(results_path)
    current_path=results_path+'subset'+str(i)
    if  not os.path.exists(current_path):
        cmd='mkdir subset'+str(i)
        os.system(cmd)
    os.chdir(results_path+'subset'+str(i))
    cmd='mkdir half1'
    os.system(cmd)
    cmd='mkdir half2'
    os.system(cmd)
    # df_1=df.sample(frac=0.5)
    transfer_fast5(total_fl[:int(total_fl.shape[0]/2)], df, results_path+'subset'+str(i)+'/half1')
    # df_2 = pd.concat([df, df_1, df_1]).drop_duplicates(keep=False)
    transfer_fast5(total_fl[int(total_fl.shape[0]/2):], df, results_path + 'subset' + str(i) + '/half2')
    print(1)
# df.columns = [1]

