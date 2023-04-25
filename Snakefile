chroms  = [f"chr{x}" for x in range(1, 23)]

wildcard_constraints:
    chrom="chr(\d+|[XY])"

rule add_annotation_info:
    input:
        "data/anno/all_PASS_{chrom}.hg38_multianno.txt",
        "data/macs/{chrom}_snps.tsv", "data/macs/{chrom}_indels.tsv"
    output: "data/done/{chrom}.tsv"
    script: "scripts/add_annotation_info.py"

rule create_sqlite_database:
    input: expand("data/done/{chrom}.tsv", chrom=chroms)
    output: "danmac5.sqlite"
    script: "scripts/create_database.py"
