# TODO rename this file to downloads.py or so...
from snakemake.remote import FTP, HTTP

FTP = FTP.RemoteProvider()
HTTP = HTTP.RemoteProvider()

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

rule gencode_gff:
    input:
        FTP.remote("ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.annotation.gff3.gz", static=True, keep_local=True)
    output:
        "ref/gencode.vM20.annotation.gff3"
    shell:
        "gzip -d -c {input} > {output}"

rule gtrnadb:
    input:
        HTTP.remote("gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz", insecure=True, static=True, keep_local=True)
    output:
        "ref/mm10-tRNAs.gff3"
    shell:
        """
        tar -xzf {input} mm10-tRNAs.bed
        awk 'BEGIN{{OFS="\t"}}{{ print $1,"gtrnadb_awk","exon",$2+1,$3,$5,$6,".","gene_type=tRNA;gene_id="$4";gene_name="$4 }}' mm10-tRNAs.bed > {output}
        rm mm10-tRNAs.bed
        """

# grep processing is used to delete duplicate entries. and to rename miRNAs to exons as required by featureCounts (to make it compliant with the other refs)
rule mirbase:
    input:
        FTP.remote("mirbase.org/pub/mirbase/21/genomes/mmu.gff3", static=True, keep_local=True)
    output:
        "ref/mirbase_mmu21.gff3"
    shell:
        """
        grep -vP 'MIMAT\d+_\d+' {input} | sed -e 's/\<miRNA\>/exon/' -e 's/\<Name\>/gene_id/' > {output}
        """

