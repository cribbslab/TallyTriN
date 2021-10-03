__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import numpy as np
import pandas as pd
import markov_clustering as mc
from src.util.graph.bfs.ConnectedComponent import connectedComponent as gbfscc
from src.sequencing.reads.similarity.distance.Hamming import hamming
from src.util.sequence.symbol.Single import single as dnasgl
from sklearn.cluster import DBSCAN as skdbscan
from sklearn.cluster import OPTICS
from sklearn.cluster import Birch


class dbscan(object):

    def __init__(self, ):
        self.gbfscc = gbfscc()
        self.dnasgl = dnasgl()
        # print(self.dnasgl.get(universal=True))
        self.dna_map = self.dnasgl.todict(nucleotides=self.dnasgl.get(universal=True), reverse=False)
        self.dna_map_rev = self.dnasgl.todict(nucleotides=self.dnasgl.get(universal=True), reverse=True)
        # print(self.dna_map)

    def dfclusters(self, connected_components, df_umi_uniq_val_cnt, umi_uniq_mapped_rev):
        # print([*connected_components.values])
        df_ccs = pd.DataFrame({'cc': [*connected_components.values()]})
        df_ccs['cc_max_id'] = df_ccs['cc'].apply(lambda x: self.sort_vals(df_umi_uniq_val_cnt, x))
        df_ccs['cc_max_id2seq'] = df_ccs['cc_max_id'].apply(lambda x: umi_uniq_mapped_rev[x])
        # print(df_ccs)
        ccc = df_ccs['cc_max_id2seq'].apply(lambda x: self.reps(x)).values.tolist()
        cccc = pd.DataFrame(ccc)
        # print(cccc[[0, 1, 2]])
        d = skdbscan(eps=2.5, min_samples=2).fit(cccc)
        # d =  Birch(threshold=1.8, n_clusters=None).fit(cccc)
        asd = np.unique(d.labels_)
        asdas = np.array(d.labels_)
        # print(d.labels_)
        # return len(asd)
        return len(asd) + len(asdas[asdas == -1]), len(asdas[asdas == -1])

    def sort_vals(self, df_umi_uniq_val_cnt, cc):
        # print(df_umi_uniq_val_cnt.loc[df_umi_uniq_val_cnt.index.isin(cc)].sort_values(ascending=False).to_dict())
        xx = [*df_umi_uniq_val_cnt.loc[df_umi_uniq_val_cnt.index.isin(cc)].sort_values(ascending=False).to_dict().keys()]
        return xx[0]

    def dusters(self, df_umis):
        # print(df_umis)
        sss = pd.DataFrame.from_dict(df_umis, orient='index', columns=['umi_uniq'])
        sss['umi_uniq'].apply(lambda x: self.reps(x))
        # print(sss)
        ccc = sss['umi_uniq'].apply(lambda x: self.reps(x)).values.tolist()
        # print(ccc)
        cccc = pd.DataFrame(ccc)
        # print(cccc[[0, 1, 2]])
        d = skdbscan(eps=1.6, min_samples=1).fit(cccc)
        # d =  Birch(threshold=1.8, n_clusters=None).fit(cccc)
        asd = np.unique(d.labels_)
        asdas = np.array(d.labels_)
        # print(d.labels_)
        return len(asd) + len(asdas[asdas == -1]), len(asdas[asdas == -1])

    def reps(self, x):
        umi_ltr2digi = [self.dna_map[i] for i in list(x)]
        ids_first_pos = [i*4 for i in range(len(x))]
        ids_to_be_one = [i+j for i,j in zip(umi_ltr2digi, ids_first_pos)]
        # print(ids_to_be_one)
        one_hot = np.zeros(len(x)*4)
        one_hot[ids_to_be_one] = 1
        return one_hot.astype(np.int)