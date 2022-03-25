dbsnp138_vcf_url = (
    "https://storage.googleapis.com/genomics-public-data/references/"
    "hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
)

rule get_dbsnp_vcf:
    output: "resources/hg38.dbsnp138.vcf"
    params: url=dbsnp138_vcf_url
    shell: "wget {params.url} -O {output}"

rule unpack_tarballs:
    input: 
        "00-tarballs/DNK_Hansen_Nyegaard_Werge_cRAM_VCF_sumstat_all.tar.gz",
        "00-tarballs/DNK_Hansen_Nyegaard_Werge_cRAM_VCF_sumstat_female.tar.gz",
        "00-tarballs/DNK_Hansen_Nyegaard_Werge_cRAM_VCF_sumstat_male.tar.gz"
    output: directory("01-unpacked/")
    shell: """
        mkdir {output}
        for f in {input}; do
            tar -C {output} -xf $f --strip-components 5 ;
        done
    """

wildcard_constraints:
    chrom="chr(\d+|[XY])"

rule process_vcfs:
    input: 
        "01-unpacked/{chrom}.all.sumstats.raremasked.vcf.gz",
        "01-unpacked/{chrom}.female.sumstats.raremasked.vcf.gz",
        "01-unpacked/{chrom}.male.sumstats.raremasked.vcf.gz"
    output: "02-processed/{chrom}.tsv"
    script: "scripts/process_vcfs.R"

chroms = [f"chr{x}" for x in ["X", "Y", *range(1, 23)]]

rule create_sqlite_database:
    input: 
        expand(rules.process_vcfs.output[0], chrom=chroms),
        rules.get_dbsnp_vcf.output
    output: "danmac5.db"
    shell: "sqlite3 {output} < scripts/prepare_db.sqlite3"
