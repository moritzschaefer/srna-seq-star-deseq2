import pandas as pd
import re
import copy

reads = pd.read_csv(snakemake.input.counts, sep="\t")
read_columns = copy.copy(reads.columns)

annotation = pd.read_csv(
    snakemake.input.annotation,
    comment='#',
    sep='\t',
    names=[
        "seq_id", "source", "type", "start", "end", "score", "strand", "phase",
        "attributes"
    ])

annotation['gene_id'] = annotation['attributes'].apply(
    lambda v: re.search('gene_id=([^;]*)', v).groups()[0])
annotation['gene_type'] = annotation['attributes'].apply(
    lambda v: re.search('gene_type=([^;]*)', v).groups()[0])
annotation['gene_name'] = annotation['attributes'].apply(
    lambda v: re.search('gene_name=([^;]*)', v).groups()[0])

reads = reads.merge(
    annotation[['gene_id', 'gene_type',
                'gene_name']].groupby('gene_id').first().reset_index(),
    left_on='Geneid',
    right_on='gene_id')

type_filters = {'miRNA': lambda v: v['Geneid'].str.contains('Mir')}

for rna_type in ['snRNA', 'snoRNA', 'Mt_tRNA']:
    sub_reads = reads[reads['gene_type'] == rna_type].copy()
    sub_reads['Geneid'] = sub_reads['gene_name']
    sub_reads['Geneid'] = sub_reads['gene_name']
    # avoid duplicates
    sub_reads = sub_reads[read_columns].groupby('Geneid').first().reset_index()
    sub_reads.to_csv(
        snakemake.output[f'{rna_type.lower()}s'], index=False, sep='\t')
