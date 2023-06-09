import numpy as np
import pandas as pd
from collections import OrderedDict
import argparse
# Created in 20230323
# Created by QUO(PH)
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
                    fa_dic[full_name] = "".join(seq_l)  # store previous one

                full_name = line.strip().replace(">", "")
                short_name = full_name.split(" ")[0]
                seq_l = []
            else:  # collect the seq lines
                if len(line) > 8:  # min for fasta file is usually larger than 8
                    seq_line1 = line.strip()
                    seq_l.append(seq_line1)

        fa_dic[short_name] = "".join(seq_l)  # store the last one
    return fa_dic


def extract_tombo(ab_path,fasta_path):
    tombo_result_minus_path= ab_path + ".level_samp_comp_detect.statistic.minus.wig"
    tombo_result_plus_path= ab_path +".level_samp_comp_detect.statistic.plus.wig"
    tombo_difference_minus_path= ab_path + ".level_samp_comp_detect.difference.minus.wig"
    tombo_difference_plus_path= ab_path + ".level_samp_comp_detect.difference.plus.wig"
    chrome_name=list(read_fasta_to_dic(fasta_path).keys())[0]
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
        df.insert(loc=0, column="A", value=np.repeat([[chrome_name]],df.shape[0]))
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
    df = pd.concat([plus_df, minus_df], ignore_index=True)
    all_df=pd.merge(df,df_new,how='inner',on=["A",0,"A1","A2",3,2])
    all_df.columns=['A', 0, 'A1', 'A2', '-log10(P_value)', 3 ,2, 'difference']
    all_df=all_df[['A', 0, 'A1', 'A2', 3,2,'-log10(P_value)','difference']]
    return all_df
def main(args):
    for i in range(1,101):
        total_path=args.input+'subset'+str(i)+'/sample'
        fasta_path=args.fasta
        try:
            df=extract_tombo(total_path,fasta_path)
            df.to_csv(total_path+"_tombo.bed",sep="\t",index=None,header=None)
            print("saved as "+total_path+"_tombo.bed")
        except Exception:
            print(Exception)


if __name__ == '__main__':
    # global npd_file
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default="SIN/subsets/",
                        help='The path and suffix of tombo result')
    parser.add_argument('--fasta', default="SINV_Toto1101.fa",
                        help='reference path')
    args = parser.parse_args()
    main(args)