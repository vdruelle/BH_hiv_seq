rule sub_sample:
    message:
        """
        Sub sampling the raw data using seqtk.
        """
    input:
        sequences = "data/raw/{region}.fasta"
    output:
        sequences = "data/raw/{region}_subsampled.fasta"
    params:
        nb_sequences = 1000
    shell:
        """
        seqtk sample -s100 {input.sequences} {params.nb_sequences} > {output.sequences}
        """

rule metadata:
    message:
        "Creating metadata file from sequences names."
    input:
        sequences = rules.sub_sample.output.sequences
    output:
        metadata = "data/raw/{region}_metadata.tsv"
    shell:
        """
        python scripts/metadata_from_names.py {input.sequences} {output.metadata}
        """


rule align:
    message:
        "Aligning sequences to {input.reference} using Augur. Strip gaps relative to reference."
    input:
        sequences = rules.sub_sample.output.sequences,
        reference = "data/reference/HXB2_{region}.fasta"
    output:
        alignment = "data/alignments/to_HXB2/{region}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
        """

rule tree:
    message:
        "Building tree using augur"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "intermediate_files/tree_{region}.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """



rule clean:
    message: "Removing generated files."
    shell:
        "rm data/raw/*subsampled.fasta"
        "rm data/alignments/to_HXB2/*"


# rule lanl_to_augur:
#     message:
#         "Generating augur format metadata from lanl metadata"
#     input:
#         metadata = "data/raw/{region}_metadata.tsv"
#     output:
#         metadata = "data/raw/{region}_metadata_augur.tsv"
#     shell:
#         "python scripts/lanl_to_augur_metadata.py {input.metadata}"
