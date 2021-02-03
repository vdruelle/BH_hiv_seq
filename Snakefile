def lanl_to_augur_metadata(lanl_file):
    "Transform the lanl metadata to augur format"

    import pandas as pd
    df = pd.read_csv(lanl_file, sep="\t", skiprows=2)
    df.head()


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
        Sub sampling the raw data using augur filter.
        """
    input:
        sequences = "data/raw/{region}.fasta",
        metadata = "data/raw/{region}_metadata_augur.tsv"
    output:
        sequences = "data/raw/{region}_subsampled.fasta"
    params:
        nb_sequences = 1000
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --sequences-per-group {params.nb_sequences} \
            --output {output.sequences}
        """


rule clean:
    message: "Removing generated files."
    shell:
        "rm data/raw/*subsampled.fasta"
