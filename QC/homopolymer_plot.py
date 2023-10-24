import pandas as pd
from plotnine import *
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
def plot_observe_homo(input_path):
    df2 = pd.read_csv(input_path, sep="\t")
    df2["value"] = df2["value"] * 100

    Homo = (
            ggplot(df2, aes(x="base", y="value", group="group",color='group')) +
            geom_line(size=1) +
            geom_point(size=3) +
            scale_colour_manual(values=["#8ECFC9",'#FFBE7A', '#FA7F6F'])+
            scale_fill_manual(values=["#8ECFC9", '#FFBE7A', '#FA7F6F']) +
            xlab("Base") + theme_bw() +
            labs(y="Accuracy of homopolymer identification (%)")+
            theme(axis_text=element_text(size=12, color="black",family='Arial'),
                  axis_title=element_text(size=12, color="black",family='Arial'))

    )
    print(Homo)
    Homo.save(filename="plots/Observed_homopolymer_identification.pdf",
              width=5, height=4, dpi=300)
plot_observe_homo("./data/homo_new_basecall.txt")