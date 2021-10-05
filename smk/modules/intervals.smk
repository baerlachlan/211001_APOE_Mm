rule intervals:
    output:
        "refs/exons.intervals"
    params:
        species = config.species,
        ensembl_release = config.ensembl_release
    conda:
        "../envs/R.yaml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-00:30:00"
    script:
        "intervals.R"