rule all:
    input:
        auspice_json = "visualisation/pol.json"

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
        metadata = "data/raw/{region}_subsampled_metadata.tsv"
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

rule split_positions:
    message:
        "Splitting {input.alignment} into 1st, 2nd and 3rd position alignments."
    input:
        alignment = rules.align.output.alignment
    output:
        alignment_first = "data/alignments/to_HXB2/{region}_1st.fasta"
        alignment_second = "data/alignments/to_HXB2/{region}_2nd.fasta"
        alignment_third = "data/alignments/to_HXB2/{region}_3rd.fasta"
    shell:
        """
        python scripts/split_positions.py {input.alignment}
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

rule refine:
    message:
        """
        Refining tree for timetree
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = rules.metadata.output.metadata
    output:
        tree = "intermediate_files/timetree_{region}.nwk",
        node_data = "intermediate_files/branch_lengths_{region}.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment
    output:
        node_data = "intermediate_files/{region}_nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule export:
    message: "Exporting data files for visualisation in auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.metadata.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
    output:
        auspice_json = "visualisation/{region}.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} \
            --output {output.auspice_json} \
            --title HIV-1_{wildcards.region}
        """


rule clean:
    message: "Removing generated files."
    shell:
        """
        rm data/raw/*subsampled.fasta
        rm data/alignments/to_HXB2/*
        rm intermediate_files/*
        rm visualisation/*
        """
