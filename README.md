# Nextflow-RNASEQ
This is a project for RNA Sequence Analysis using Nextflow pipeline.

# Environment Setup
- install nextflow
    ```bash
    cd ~
    wget -qO- https://get.nextflow.io | bash
    chmod +x nextflow
    ```
    - You can move the executable nextflow to /usr/bin/ to directly use it
- Dependencies
    - conda
    - TO BE TEST ON ERISONE

# Pipeline Confioguratio
- For basic configuration, reference to ***./conf/sample.config*** or ***./conf/test.config***
    - input/output path
    - genome type
        - preconfigured genomes saved in ***./conf/igenomes.config***
        - you can choose genome type from ***./conf/igenomes.config***
        - you can also provide ***STAR Index***, ***gtf*** by your self

- For advanced configuration, reference to ***./nextflow.config***
    - like parameters for each algorithm or tool
    - check the parameters in the file
    - config them in ***./conf/YOUR_SAMPLE_NAME.config***
    - you can reference to ***./conf/sample.config*** -> **STAR Options**
