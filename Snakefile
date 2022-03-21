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

dbsnp138_vcf_url = (
    "https://storage.googleapis.com/genomics-public-data/references/"
    "hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
)

rule get_dbsnp_vcf:
    output: "resources/hg38.dbsnp138.vcf"
    params: url=dbsnp138_vcf_url
    shell: "wget {params.url} -O {output}"

rule create_sqlite_database:
    input: 
        rules.unpack_tarballs.output, 
        rules.get_dbsnp_vcf.output
    output: "danmac5.db"
    shell: "sqlite3 {output} < scripts/prepare_db.sqlite3"
