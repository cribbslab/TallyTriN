from src.util.sequence.convert.ManyToSingle import manyToSingle as sm2s
from src.sequencing.bulk.reads.simulate.Native import native as simunat

# download cdna lib at http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# define the data path
DEFINE = {
    'fasta_fpn': 'data/omics/genomics/fasta/cdna/GRCh38/Homo_sapiens.GRCh38.cdna.all.fa',
    'sv_fpn_ids': 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n.txt',
    'sv_fasta_fp': 'data/omics/genomics/fasta/cdna/GRCh38/',
}

# ## generate ENST map
sm2s(
    fasta_fpn=DEFINE['fasta_fpn']
).svid(
    sv_fpn=DEFINE['sv_fpn_ids']
)

# ## split the whole fasta into individual seqs and optimize I/O
sm2s(
    fasta_fpn=DEFINE['fasta_fpn']
).save(
    sv_fp=DEFINE['sv_fasta_fp']
)

# ## simulate sequences by the native method
# this will create a seq arr consisting of cdna_num (e.g., 1000) sequences
simu_by_native = simunat(
    cand_pool_fpn=DEFINE['sv_fpn_ids'],
    cdna_fp=DEFINE['sv_fasta_fp'],
    cdna_num=1000
)
seqs = simu_by_native.generate()