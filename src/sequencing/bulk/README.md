# simulate

## define the data path
```html
DEFINE = {
    'fasta_fpn': 'data/omics/genomics/fasta/cdna/GRCh38/Homo_sapiens.GRCh38.cdna.all.fa',
    'sv_fpn_ids': 'data/omics/genomics/fasta/cdna/GRCh38/cdna_n.txt',
    'sv_fasta_fp': 'data/omics/genomics/fasta/cdna/GRCh38/',
}
```

## generate ENST map
```html
sm2s(
    fasta_fpn=DEFINE['fasta_fpn']
).svid(
    sv_fpn=DEFINE['sv_fpn_ids']
)
```

## split the whole fasta into individual seqs and optimize I/O
```html
sm2s(
    fasta_fpn=DEFINE['fasta_fpn']
).save(
    sv_fp=DEFINE['sv_fasta_fp']
)
```

## simulate sequences by the native method - this will create a seq arr consisting of cdna_num (e.g., 1000) sequences
```html
simu_by_native = simunat(
    cand_pool_fpn=DEFINE['sv_fpn_ids'],
    cdna_fp=DEFINE['sv_fasta_fp'],
    cdna_num=1000
)
seqs = simu_by_native.generate()

```