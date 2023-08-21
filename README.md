¿Prefieres el español? Puedes encontrar el README en español [aquí](README.es.md).

# nCoV-2019 Analysis Pipeline (Oxford Nanopore)

Nextflow pipeline for assembly and underlying analysis of SARS-CoV-2 genomes.
Uses the ARTIC bioinformatic protocol, so using their sequencing protocol is required ([more info](https://artic.network/ncov-2019)).

Developed to automate the bioinformatic analysis associated to the virus genomic surveillance initiative pushed by the chilean government and executed by various teams in the country, like us at the University of Magallanes, in the city of Punta Arenas. More information about us in [catg.cl](https://catg.cl) (link in spanish).

Currently includes:

- Assembly and variant calling using [ARTIC](https://github.com/artic-network/fieldbioinformatics)
- Variant annotation using [SnpEff](https://pcingola.github.io/SnpEff/)
- Identification of clades/lineages of Pangolin and Nexstrain (including common variant name)
- Quality metrics and statistics (qc, coverage)
- Generation of final summary in Excel format containing relevant information obtained by the pipeline
- Inicial preparation for GISAID submission

## Quick Start

1. Install Nextflow and Apptainer (or Docker if you prefer that).
2. Create a CSV file with the columns `sample` and `fastq_file`, for this a single FASTQ file is needed for each sample. If you plan to submit your results to GISAID, it is recommended to have a third column: `gisaid_covv_virus_name` with the names to use for each sample in the submission (the format must be `hCoV-19/<Country>/<ID>/<Year>`).
3. Execute the pipeline. Examples:
   Apptainer:

   ```
   nextflow run catg-umag/ncov2019-ont-nf -r v2.0.0 -profile apptainer \
       --sample_data <csv_file>
   ```

   Docker, with nanopolish instead of medaka:

   ```
   nextflow run catg-umag/ncov2019-ont-nf -r v2.0.0 -profile docker \
       --sample_data <csv_file> --artic_use_nanopolish \
       --fast5_directory <fast5_dir> --sequencing_summary <txt_file>
   ```

4. When the execution finishes, you will find your results in the `results/` directory.

## Requirements

- [Nextflow](https://www.nextflow.io/) (>= 22.04.0)
- [Apptainer](https://apptainer.org/docs/admin/main/installation.html) or [Docker](https://www.docker.com/get-started)

Nextflow works best in GNU/Linux and MacOS. If you really must use Windows, you need to setup everything on [WSL](https://learn.microsoft.com/en-us/windows/wsl/install).

## Pipeline Utilization

### Inputs

This pipeline is designed to run from the data already basecalled and demultiplexed, so the FASTQ files. If you don't have this yet or want to make basecalling again with higher accuracy, check [this pipeline](https://github.com/catg-umag/ont-basecalling-demultiplexing).

So, you will need:
- a CSV file the information for each sample (more on that later)
- a single FASTQ file for each sample (compressed or not). If you have the default files generated with Guppy (lots of files), you should concatenate them first.
- subdirectory `pass` of the directory with the FAST5 files, and `sequencing_summary.txt` (only if you want to use Nanopolish instead of Medaka)


### How to prepare your sample information file

This file (in CSV format) must have at least two columns: `sample` and `fastq_file`, associating each sample with their corresponding FASTQ data. The values of the sample columns will be used to name files associated to each sample, so it's recommended to keep it simple (ideally only letters and/or numbers). You can also add any additional column you want, these will be included in the final summary. It is recommended to use absolute paths for the FASTQ files, but you can also use relative paths (relative to the execution directory).

Another set of columns you can include in this file are the related to the GISAID submission. These can be any of the columns in the submission template (available in `data/`), indicated by the `gisaid_` prefix followed by the column code (the first row of the Excel file). It is recommended to at least include `covv_virus_name`, because the pipeline will use this column to rename the sequences in the FASTA file, process that you will have to do manually if you don't provide these values. You can find in the Excel file all the available fields and their required formats. Some columns will be filled by the pipeline either way, these are: `fn`, `covv_type`, `covv_passage`, `covv_host`, `covv_seq_technology`, `covv_assembly_method` and `covv_coverage`.

For example for the following 5 samples we included the extra field `city` (that will be available in the final summary), and the GISAID fields `covv_virus_name`, `covv_collection_date` and `covv_location`.

```
sample,fastq_file,city,gisaid_covv_virus_name,gisaid_covv_collection_date,gisaid_covv_location
2101011,./data_fastq/barcode01.fastq.gz,Punta Arenas,hCoV-19/Chile/CADIUMAG-51/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2101018,./data_fastq/barcode02.fastq.gz,Punta Arenas,hCoV-19/Chile/CADIUMAG-52/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2101053,./data_fastq/barcode03.fastq.gz,Punta Arenas,hCoV-19/Chile/CADIUMAG-53/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2101025,./data_fastq/barcode04.fastq.gz,Punta Arenas,hCoV-19/Chile/CADIUMAG-54/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2105001,./data_fastq/barcode05.fastq.gz,Punta Arenas,hCoV-19/Chile/CADIUMAG-55/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
```

### Pipeline parameters

The pipeline has various parameters to help to suit your needs. These are:

| Parameter                      | Required | Default                      | Description                                                                                                                                             |
| ------------------------------ | -------- | ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| sample_data                    | yes      | ---                          | CSV file with each sample related information.                                                                                                          |
| fast5_directory                | no       | ---                          | FAST5 `pass` subdirectory, only needed if you want to use Nanopolish.                                                                                   |
| sequencing_summary             | no       | ---                          | Summary of the sequencing process, required if you provide `fast5_directory`.                                                                           |
| run_id                         | no       | ""                           | Optional ID that will be used as suffix in the summary files.                                                                                           |
| artic_primer_scheme            | no       | SARS-CoV-2/V4.1              | ARTIC primer scheme used.                                                                                                                               |
| artic_medaka_model             | no       | r103_hac_g507                | Medaka model to use (this needs to match the basecalling model used).                                                                                   |
| artic_use_nanopolish           | no       | false                        | Use Nanopolish instead of Mekada for polishing.                                                                                                         |
| artic_normalise                | no       | 500                          | Coverage target used by the ARTIC pipeline.                                                                                                             |
| gisaid_template                | no       | data/20230515_EpiCoV....xlsx | Path to the GISAID upload template (included in the repository).                                                                                        |
| gisaid_submission_enabled      | no       | true                         | Enable (or not) the GISAID submission preparation.                                                                                                      |
| publish_minimum_completion     | no       | 95                           | Percentage value (0 - 100) indicating the required percentage of covered bases of the reference to call an assembly "valid" for the submission process. |
| multiqc_config                 | no       | conf/multiqc_config.yaml     | Path to a valid MultiQC config to use. The pipeline includes one already                                                                                |
| output_directory               | no       | results                      | Directory to storage the pipeline results.                                                                                                              |

These parameters / default values are defined in `params.default.yml`.
You can specify these by command line options (for example `--run_id 20210501A`), but you can also provide them in a YAML file, use the `params.default.yml` as template (don't edit or move the file, is required for the execution).

### Execution

The pipeline can be downloaded directly (with `git clone` for example), and you can specify to Nextflow the path, but you can also indicate just `catg-umag/ncov2019-ont-nf` and a version with the `-r` option, and that will download the pipeline automatically. For example:

```
nextflow run catg-umag/ncov2019-ont-nf -r v2.0.0 -profile apptainer ...
```

Is mandatory to use either the `apptainer` or `docker` profiles, because if you don't the pipeline will expect to have all the necessary tools installed and most probably won't work.

Example using parameters by command line:

```
nextflow run ncov2019-ont-nf/ -profile apptainer --sample_data input/samples.csv --run_id R210505
```

Example using a YAML file for the parameters:

`params.yml`:

```yaml
sample_data: input/samples.csv
run_id: R210505
artic_medaka_model: r941_min_hac_g507
artic_normalise: 1000
gisaid_submission_enabled: false
```

```
nextflow run ncov2019-ont-nf/ -profile apptainer -params-file params.yml
```

By default the pipeline will execute locally, but if you are using a Slurm Cluster, you can use the `slurm` profile. And don't forget the `apptainer` profile too (this must be installed in the cluster as well): `-profile slum,apptainer`.

## Results

Inside the results directory, you will find the following:

- `raw_data/`: FASTQ files for each sample, generated after basecalling+demultiplexing and length filter
- `artic/`: results of the ARTIC pipeline execution for each of the samples
- `qc/`: quality and coverage metrics for each of the samples
- `vcf/`: VCF files annotated by SnpEff
- `lineages/`: RAW results from the lineage identification tools
- `summary/`: directory with summaries of the obtained information:
  - `all_consensus.fasta`: FASTA file with the consensus for each sample (even the incomplete ones)
  - `sample_summary.csv`: CSV file with the provided info for each sample and also lineages and coverage stats
  - `variants_list.csv`: CSV with the variants for all the samples with their relevant info
  - `summary.xlsx`: Excel file with the information of `sample_summary.csv` as first sheet and a summary of the variant presence on the samples in the second one.
- `gisaid_submission/`: if you enabled the submission preparation you will find here a FASTA file filtered by coverage and with sequences renamed if you provided the corresponding name, and also the upload template filled with the available information (you could still need to fill some values manually)

## How to Cite

If this work was useful to you, you can cite this article:

> González-Puelma, J.; Aldridge, J.; Montes de Oca, M.; Pinto, M.; Uribe-Paredes, R.; Fernández-Goycoolea, J.; Alvarez-Saravia, D.; Álvarez, H.; Encina, G.; Weitzel, T.; Muñoz, R.; Olivera-Nappa, Á.; Pantano, S.; Navarrete, M.A. Mutation in a SARS-CoV-2 Haplotype from Sub-Antarctic Chile Reveals New Insights into the Spike’s Dynamics. Viruses 2021, 13, 883. https://doi.org/10.3390/v13050883

## Acknowledgements

- ARTIC and all its dependencies
- Python and libraries: biopython, openpyxl, pandas, scipy
- Bioinformatic tools: FastQC, MultiQC, samtools, SnpEff, Nextclade, Pangolin
- Other tools: Nextflow, libreoffice
