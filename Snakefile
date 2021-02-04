rule lanl_to_augur:
    message:
        "Generating augur format metadata from lanl metadata"
    input:
        metadata = "data/raw/{region}_metadata.tsv"
    output:
        metadata = "data/raw/{region}_metadata_augur.tsv"
    shell:
        "python scripts/lanl_to_augur_metadata.py {input.metadata}"


rule sub_sample:
    message:
        """
        Sub sampling the raw data using seqtk. Adding HXB2 sequence.
        """
    input:
        sequences = "data/raw/{region}.fasta",
        HXB2 = "data/reference/HXB2.fasta"
    output:
        sequences = "data/raw/{region}_subsampled.fasta"
    params:
        nb_sequences = 1000
    shell:
        """
        seqtk sample -s100 {input.sequences} {params.nb_sequences} > {output.sequences}
        python scripts/add_HXB2.py {output.sequences} {input.HXB2}
        """



rule clean:
    message: "Removing generated files."
    shell:
        "rm data/raw/*subsampled.fasta"
