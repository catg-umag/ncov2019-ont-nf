¿Prefieres el español? Puedes encontrar el README en español [aquí](README.es.md).

# nCoV-2019 Analysis Pipeline (Oxford Nanopore Sequencing)

Nextflow pipeline for assembly and underlying analysis of SARS-CoV-2 genomes.
Uses ARTIC's bioinformatic protocol, so using its protocol for the sequencing step is required [more info](https://artic.network/ncov-2019).

Created to automate the bioinformatic analysis associated to the virus genomic surveillance initiative pushed by the chilean government and executed by various teams in the country, like us, at the University of Magallanes, in the city of Punta Arenas.

Currently includes:

- Basecalling and demultiplexing using ONT Guppy (optional, requires Guppy previously installed)
- Assembly and variant calling using [ARTIC](https://github.com/artic-network/fieldbioinformatics)
- Variant annotation using [SnpEff](https://pcingola.github.io/SnpEff/)
- Identification of clades/lineages of GISAID, Pangolin and Nexstrain
- Quality metrics and statistics (qc, coverage)
- Generating of final summary in Excel format containing relevant information obtained by the pipeline
- Inicial preparation for GISAID submission

## Quick Start

1. Install Nextflow, Singularity (Docker is an option too) and Guppy (only if basecalling/demux needed)
2. Create a CSV file the columns `barcode` and `sample`, specifying a sample name for each barcode used. If you plan to submit your results to GISAID, it is recommended to have a third column: `gisaid_covv_virus_name` having the names to use for each sample in the submission (the format must be `hCoV-19/<Country>/<ID>/<Year>`).
3. Execute the pipeline. Examples:
   With basecalling, singularity:

   ```
   nextflow run catg-umag/ncov2019-ont-nf -r v1.0 -profile singularity \
       --sample_data <csv_file> --fast5_directory <fast5_dir>
   ```

   No basecalling, docker

   ```
   nextflow run catg-umag/ncov2019-ont-nf -r v1.0 -profile docker \
       --sample_data <csv_file> --fast5_directory <fast5_dir> \
       --fastq_directory <fastq_dir> --sequencing_summary <txt_file>
   ```

4. When the execution finishes, you will find your results in the `results/` directory.

# Requirements

- [Nextflow](https://www.nextflow.io/) (>= 20.07)
- [Singularity](https://sylabs.io/guides/3.7/admin-guide/) o [Docker](https://www.docker.com/get-started)
- Guppy (only if basecalling needed, available in the [Oxford Nanopore Community](https://community.nanoporetech.com/downloads))

Nextflow works best in GNU/Linux and MacOS. If you really want to use Windows, you need to setup everything on [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10)

## Pipeline Utilization

### Inputs

Depending your previous steps we have two options:

1. You already did basecalling and demultiplexing (for example, with MinKNOW while sequencing), and you don't want a new one. In this case, you will need:
   - subdirectory `pass` of the directory with the FAST5 files
   - directory with the data basecalled and demultiplexed, this directory should have the `barcodeXX` subdirectories
   - the `sequencing_summary.txt` file, you should find it somewhere in the FAST5 and/or FASTQ directory
   - a CSV file the information for each sample (more on that later)
2. You want basecalling as part of the pipeline, because yo didn't do it previously or because you want the high accuracy one. Consider that the high accuracy basecalling is really expensive computationally, and you would most probably need a supported GPU. In this case, you need:
   - the whole `fast5` directory
   - the CSV file with the samples information

### How to prepare your sample information file

This file (in CSV format) must have at least two columns: `barcode` and `sample`, associating each barcode to a sample. The values of the sample columns will be used to name files associated to each sample, so it's recommended to keep it simple (ideally only letters and/or numbers). You can also add every additional column you want, these will be included in the final summary.

Another set of columns you can include in this file are the related to the GISAID submission. These can be any of the columns in the submission template (available in `data/`), indicated by the `gisaid_` prefix followed by the column code (the first row of the Excel file). It is recommended to at least include `covv_virus_name`, because the pipeline will use this column to rename the sequences in the FASTA file, process that you will have to do it manually if you don't provide these values. You can find in the Excel file all the available field and their required formats. Some columns will be filled by the pipeline either way, these are: `fn`, `covv_type`, `covv_passage`, `covv_host`, `covv_seq_technology`, `covv_assembly_method` and `covv_coverage`.

For example for the following 5 samples we included the extra field `city` (that will be available in the final summary), and the GISAID fields `covv_virus_name`, `covv_collection_date` and `covv_location`.

```
sample,barcode,city,gisaid_covv_virus_name,gisaid_covv_collection_date,gisaid_covv_location
2101011,barcode01,Punta Arenas,hCoV-19/Chile/CADIUMAG-51/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2101018,barcode02,Punta Arenas,hCoV-19/Chile/CADIUMAG-52/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2101053,barcode03,Punta Arenas,hCoV-19/Chile/CADIUMAG-53/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2101025,barcode04,Punta Arenas,hCoV-19/Chile/CADIUMAG-54/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2105001,barcode05,Punta Arenas,hCoV-19/Chile/CADIUMAG-55/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
```

### Pipeline Parameters

The pipeline has various parameters to help to suit your needs. These are:

| Parameter                  | Required | Default                      | Description                                                                                                                                             |
| -------------------------- | -------- | ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| sample_data               | yes      | ---                          | CSV file with each sample related information.                                                                                                          |
| fast5_directory            | yes      | ---                          | FAST5 directory, if you don't want basecalling, point only the `pass/` subdirectory.                                                                    |
| fastq_directory            | no       | ---                          | FASTQ directory the the basecalled/demultiplexed files if you don't want a new basecalling.                                                             |
| sequencing_summary         | no       | ---                          | Summary of the sequencing process, required if you provide `fastq_directory`.                                                                           |
| run_id                     | no       | ""                           | Optional ID that will be used as suffix in the summary files.                                                                                           |
| guppy_basecalling_config   | no       | dna_r9.4.1_450bps_hac.cfg    | Basecalling configuration to use in Guppy, more on that in the Guppy manual.                                                                            |
| guppy_barcodes             | no       | "EXP-NBD104 EXP-NBD114"      | Utilized barcoding kit(s) in the sequencing process.                                                                                                    |
| guppy_demu_both_ends       | no       | true                         | Require both barcodes (5' and 3') for the demultiplexing.                                                                                               |
| guppy_enable_gpu           | no       | true                         | Use GPU for the basecalling.                                                                                                                            |
| guppy_gpu_config           | no       | "--device auto --gpu..."     | GPU basecalling configuration. The default value works well with a pair Tesla V100, but if you have something different, you could need som changes.    |
| artic_primers_scheme       | no       | nCoV-2019/V3                 | ARTIC primer scheme used.                                                                                                                               |
| artic_normalise            | no       | 500                          | Coverage target used by the ARTIC pipeline.                                                                                                             |
| gisaid_clades              | no       | data/gisaid_clades.csv       | Path of the GISAID clades specification file (included in the repository).                                                                              |
| gisaid_template            | no       | data/20210222_EpiCoV....xlsx | Path to the GISAID upload template (included in the repository).                                                                                        |
| gisaid_submission_enabled  | no       | true                         | Enable (or not) the GISAID submission preparation.                                                                                                      |
| publish_minimum_completion | no       | 95                           | Percentage value (0 - 100) indicating the required percentage of covered bases of the reference to call an assembly "valid" for the submission process. |
| output_directory           | no       | results                      | Directory to storage the pipeline results.                                                                                                              |

These parámeters / default values are defined in `params.default.yml`.
You can specify these by command line options (for example `--run_id 20210501A`), but you can also provide them in a YAML file, use the `params.default.yml` as template (don't edit or move the file, is required for the execution).

### Execution

The pipeline can be downloaded directly (with `git clone` for example), and you can specify to Nextflow the path, but you can also indicate just `catg-umag/ncov2019-ont-nf` and a version with the `-r` option, and that will download the pipeline automatically. For example:

```
nextflow run catg-umag/ncov2019-ont-nf -r v1.0 -profile singularity ...
```

Is mandatory to use either the `singularity` or `docker` profiles, because if you don't it will expect to have all the tools installed and most probably won't work.

Example using parameters by command line:

```
nextflow run ncov2019-ont-nf/ -profile singularity --fast5_directory input/fast5/ \
    --sample_data input/samples.csv --run_id R210505
```

Example using a YAML file for the parameters:

`params.yml`:

```yaml
sample_data: input/samples.csv
fast5_directory: input/fast5/
fastq_directory: input/demultiplexed/
sequencing_summary: input/sequencing_summary.txt
artic_normalise: 1000
gisaid_submission_enabled: false
```

```
nextflow run ncov2019-ont-nf/ -profile singularity -params-file params.yml
```

By default the pipeline will execute locally, but if you are using a Slurm Cluster, you can use the `slurm` profile. Ant don't forget the `singularity` profile too (this must be installed in the cluster as well): `-profile slum,singularity`.

## Results

Inside the results directory, you will find the following:

- `raw_data`: FASTQ files for each sample, generated after basecalling+demultiplexing and length filter
- `artic/`: results of the ARTIC pipeline execution for each of the samples
- `qc/`: quality and coverage metrics for each of the samples
- `vcf/`: VCF files annotated by SnpEff
- `lineages/`: RAW results from the lineage identification tools
- `summary/`: directory with summaries of the obtained information:
  - `all_consensus.fasta`: FASTA file with the consensus for each sample (even the incomplete ones)
  - `sample_summary.csv`: CSV file with the provided info for each sample and also lineages and coverage stats
  - `variants_list.csv`: CSV with the variants for all the samples with their relevant info
  - `summary.xlsx`: Excel file with the information of `sample_summary.csv` as first sheet and a summary of the variant presence on the samples in the second one.
- `gisaid_submission`: if you enabled the submission preparation you will find here a FASTA file filtered by coverage and with sequences renamed if you provided the corresponding name, and also the upload template filled with the available information (you could still need to fill some values manually)

## How to cite

If this work was useful to you, you can cite this article:

> González-Puelma, J.; Aldridge, J.; Montes de Oca, M.; Pinto, M.; Uribe-Paredes, R.; Fernández-Goycoolea, J.; Alvarez-Saravia, D.; Álvarez, H.; Encina, G.; Weitzel, T.; Muñoz, R.; Olivera-Nappa, Á.; Pantano, S.; Navarrete, M.A. Mutation in a SARS-CoV-2 Haplotype from Sub-Antarctic Chile Reveals New Insights into the Spike’s Dynamics. Viruses 2021, 13, 883. https://doi.org/10.3390/v13050883

## Acknowledgements

- ARTIC and all its dependencies
- Python and libraries: biopython, openpyxl, pandas, scipy
- Bioinformatic tools: FastQC, MultiQC, samtools, SnpEff, Nextclade, Pangolin
- Other tools: Nextflow, unoconv
