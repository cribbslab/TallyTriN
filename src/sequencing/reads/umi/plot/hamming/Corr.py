import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams


class ed(object):

    def __init__(self, ):
        pass

    def fa(self, umi_pool):
        eded = []
        for i, u1 in enumerate(umi_pool[:50]):
            for j, u2 in enumerate(umi_pool[:50]):
                eded.append([i, j, hamming().general(u1, u2)])

        tr = pd.DataFrame(eded, columns=['UMI1', 'UMI2', 'Edit distance'])

        ax = sns.kdeplot(data=tr.loc[tr['Edit distance'] > 0], x='Edit distance', linewidth=4, fill=True, common_norm=False, palette="crest", alpha=0.4,)
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Tahoma']
        ax.set_xlabel('Edit distance', fontsize=13)
        ax.set_ylabel('Density', fontsize=13)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        plt.show()

        ty = tr['Edit distance'].max()
        tr['Corr'] = tr['Edit distance'].apply(lambda x: 1-x/ty)
        print(tr)

        # print(tr.loc[:100, :])
        g = sns.relplot(
            data=tr,
            x='UMI1', y='UMI2', hue='Corr', size='Corr',
            palette="vlag", hue_norm=(-1, 1), edgecolor=".7",
            height=10, sizes=(50/1.5, 250/1.5), size_norm=(-.2, .8),
        )
        # Tweak the figure to finalize
        g.set(xlabel="", ylabel="", aspect="equal")
        g.despine(left=True, bottom=True)
        g.ax.margins(.02)
        for label in g.ax.get_xticklabels():
            label.set_rotation(90)
        for artist in g.legend.legendHandles:
            artist.set_edgecolor(".7")
        plt.show()
