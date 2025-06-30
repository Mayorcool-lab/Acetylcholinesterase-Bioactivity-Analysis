configfile: "config.yaml"

rule all:
    input:
        config["merged_data"],
        config["model_output"]

rule merge_descriptors:
    input:
        fingerprints=config["fingerprints"],
        raw_data=config["raw_data"]
    output:
        merged=config["merged_data"]
    script:
        "scripts/merge_descriptors.py"

rule train_model:
    input:
        config["merged_data"]
    output:
        config["model_output"]
    script:
        "scripts/train_model.py"
rule evaluate_model:
    input:
        merged=config["merged_data"],
        model=config["model_output"]
    output:
        config["plot_output"]
    script:
        "scripts/evaluate_model.py"
