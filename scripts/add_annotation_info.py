from duckdb import sql, read_csv
import duckdb

## get annotation file
anno_raw = read_csv(snakemake.input[0], delimiter="\t", na_values = "-")
anno = sql(
    """
    select
        Chr || ':' || Start || '_' ||
            coalesce(Ref, 0) ||'_'|| coalesce(Alt, 0) AS VARID,
        regexp_extract("Gene.refGeneWithVer", '[^;]+$') as GENE,
        CLNDISDB, CLNALLELEID, CLNDN, CLNSIG,
        "Func.refGeneWithVer" as GENE_FUNCTION,
        "ExonicFunc.refGeneWithVer" as EXON_FUNCTION,
        "AAChange.refGeneWithVer" as AA_CHANGE,
        avsnp147 as RSID
    from anno_raw
    """,
    alias = "anno"
)

## get minor alelle counts
indels = read_csv(snakemake.input[1], delimiter="\t", na_values = "NA")
snps = read_csv(snakemake.input[2], delimiter="\t", na_values = "NA")
counts = sql(
    """
    select
        "#CHROM" || ':' || POS || '_' || REF || '_' || ALT  AS VARID,
        "#CHROM" AS CHROM, POS,
        "ALL"  as MAC_ALL,
        FEMALE as MAC_FEMALE,
        MALE   AS MAC_MALE
    from (
        SELECT * FROM snps
        UNION ALL
        SELECT * FROM indels
    )
    """,
    alias = "counts"
)

## merge files and write to disk
anno.join(counts, "VARID").to_csv(
    snakemake.output[0], sep="\t", header=True, na_rep=".", quotechar=""
)
