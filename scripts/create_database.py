import sqlite3
import subprocess

with sqlite3.connect(snakemake.output[0]) as con:
    con.execute(
        """
        create table danmac5 (
            VARID varchar,
            GENE varchar,
            CLNDISDB varchar,
            CLNALLELEID varchar,
            CLNDN varchar,
            CLNSIG varchar,
            GENE_FUNCTION varchar,
            EXON_FUNCTION varchar,
            AA_CHANGE varchar,
            RSID varchar,
            CHROM varchar,
            POS int,
            MAC_ALL varchar,
            MAC_FEMALE varchar,
            MAC_MALE varchar
        );
        """
    )

## BULK IMPORT ################################################################

# sqlite 'dot' commands for fast import
sqlite_commands = [".mode tabs", ".headers on", ".import '| cat -' danmac5" ]

for file in snakemake.input:
    with open(file) as f:
        proc = subprocess.run(
            ["sqlite3", snakemake.output[0], *sqlite_commands],
            stdin=f
        )

## ADD INDEXING ###############################################################

with sqlite3.connect(snakemake.output[0]) as con:
    con.executescript(
        """
        create index idx_varid on danmac5 (varid);
        create index idx_chrom_pos on danmac5 (chrom, pos);
        create index idx_chrom on danmac5 (chrom);
        create index idx_gene on danmac5 (gene);
        """
    )
