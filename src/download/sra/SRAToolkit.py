__version__ = "v1.0"
__copyright__ = "Copyright 2021"
__license__ = "MIT"
__lab__ = "Adam Cribbs lab"

import subprocess
from pysradb.sraweb import SRAweb


class sraToolkit(object):

    def __init__(self, sra_acc, todo_list):
        self.todo_list = todo_list
        self.sra_acc = sra_acc
        self.sra_db = SRAweb()
        self.df = self.sra_db.sra_metadata(self.sra_acc, detailed=True)
        self.df = self.df.loc[self.df['run_accession'].isin(self.todo_list)]
        print(self.df)

    def serial(self, ):
        # this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
        for sra_id in self.todo_list:
            print("Currently downloading: " + sra_id)
            prefetch = "prefetch " + sra_id
            print("The command used was: " + prefetch)
            subprocess.call(prefetch, shell=True)

        # this will extract the .sra files from above into a folder named 'fastq'
        for sra_id in self.todo_list:
            print("Generating fastq for: " + sra_id)
            fastq_dump = "fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/public/sra/" + sra_id + ".sra"
            print("The command used was: " + fastq_dump)
            subprocess.call(fastq_dump, shell=True)
        return 0

    def paralell(self, ):
        pass


if __name__ == "__main__":
    p = sraToolkit(
        # sra_acc='PRJNA232531',
        # todo_list=['SRR1058003', 'SRR1058023'],

        sra_acc='PRJNA274274',
        todo_list=['SRR1784310', 'SRR1784313', 'SRR1784314', 'SRR1784315'],
    )

    print(p.serial())