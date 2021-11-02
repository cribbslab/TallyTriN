__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

from pysradb.sraweb import SRAweb
from joblib import Parallel, delayed


class pysradb(object):

    def __init__(self, sra_acc, todo_list):
        self.todo_list = todo_list
        self.sra_acc = sra_acc
        self.sra_db = SRAweb()
        self.df = self.sra_db.sra_metadata(self.sra_acc, detailed=True)
        self.df = self.df.loc[self.df['run_accession'].isin(self.todo_list)]
        print(self.df.columns)

    def serial(self, ):
        return self.sra_db.download(df=self.df, url_col="sra_url_alt")

    def biparalell(self, ):
        df1, df2 = self.df.iloc[:int(self.df.shape[0]/2), :], self.df.iloc[int(self.df.shape[0]/2):, :]
        print(df1)
        print(df2)
        def single_download(df_single):
            self.sra_db.download(df=df_single, skip_confirmation=True)
        Parallel(n_jobs=2)(delayed(single_download)(df_x) for df_x in [df1, df2])
        return 0


if __name__ == "__main__":
    p = pysradb(
        # sra_acc='PRJNA232531',
        # todo_list=['SRR1058003', 'SRR1058023'],

        sra_acc='PRJNA274274',
        todo_list=['SRR1784310', 'SRR1784313', 'SRR1784314', 'SRR1784315'],
    )

    # print(p.serial())

    print(p.biparalell())