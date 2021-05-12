# Artic Sarscov2

## ¿Qué incluye?

## ¿Cómo usarlo?
### Requerimientos

Primero, se necesita tener instalado Nextflow (>=20.07) y Singularity.
En caso de querer realizar el basecalling de alta presición, necesitará contar con una GPU.

#### Preparación de los archivos de entrada
#### Ejecución del flujo de trabajo 

Para correr el pipeline ejecutar: 

` nextflow run Nanopore_SARS-CoV-2/main.nf -params-file  <params> -profile <profile>`

En el archivo <params>, se deben proporcionar los archivos de entrada, configuraciones y otras opciones. Estas son:

| Parámetro | Requerido | Por defecto | Descripción |
| ------ | ------ | --- | ---- |
| run_id | yes | --- | Identificador de la secuenciación |
| fast5 | yes | --- | Ruta al directorio FAST5 obtenido en la secuenciación. |
| fastq | no | --- | Ruta al directorio FASTQ obtenido en la secuenciación. Solo en caso de querer saltar el basecalling de alta precisión y comenzar el pipeline con los archivos fastq demultiplexados. |
| barcodes | yes | --- | Kit utilizado en la secuenciación para multiplexar las muestras. |
| sequencing-summary | yes | --- | Resumen de secuenciación obtenido en la secuenciación|
| primers_scheme | yes | nCoV-2019/V3 | Ruta al directorio con el esquema de primers utilizado en la preparación de librería. Debe estar en el directorio actual, no se recomienda usar rutas absolutas. |
| guppy_config | yes | dna_r9.4.1_450bps_hac.cfg | Configuración a utilizar en el basecalling de alta precisión. Más información sobre las configuraciones disponibles en el manual de Guppy. |
| guppy_device | no | 'cuda:0,1' | Dispositivo GPU a utilizar en el basecalling de alta precisión. |
| primers | yes | artic-ncov2019/primer_schemes | Ruta al directorio que contiene todos los esquemas de primers. Esquemas de primers disponibles en el repositorio artic-ncov2019.  |
| hac_basecalling | yes | True | Valor booleano para la ejecución del basecalling de alta precisión. False para comenzar con los archivos FASTQ demultiplexados, True para empezar con los archivos Fast5 y el basecalling de alta precisión. |
| len_min_amplicon | yes | 400 | Tamaño mínimo de los amplicones.  |
| len_max_amplicon | yes | 700 | Tamaño máximo de los amplicones más 200.  |
| ourdir | yes | results | Directorio para almacenar los resultados del pipeline. |

/*Falta samples que hay que reestructurarlo*/


Por ejemplo, una ejecución que incluya el basecalling de alta precisión y utilizando los parámtros por defecto, sería algo así:

`nextflow run Nanopore_SARS-CoV-2/main.nf --run_id 20210501 --fast5 nanopore_sequencing/fast5/ --barcodes EXP-NBD104 --sequencing-summary nanopore_sequencing/sequencing_summary_FAO03069_83da979a.txt `

Y una ejecución sin basecalling de alta precisión (comenzando con los archivos fastq demultiplexados), luciria así:

`nextflow run Nanopore_SARS-CoV-2/main.nf --run_id 20210501 --fast5 nanopore_sequencing/fast5/ --fastq nanopore_sequencing/fastq_pass/barcode0[1-6] --barcodes EXP-NBD104 --sequencing-summary nanopore_sequencing/sequencing_summary_FAO03069_83da979a.txt --hac_basecalling False`

También se puede proporcionar un archivo yaml que contenga todos los parametros a configurar, de esta manera no es necesario escribir todo en la línea de comando. Para esto descargue el archivo de ejemplo params.example.yaml y editelo de acuerdo a sus configuraciones. Luego puede ejecutar el pipeline de la siguiente manera: 

` nextflow run Nanopore_SARS-CoV-2/main.nf -params-file params.yml `


## What is included?
## How to use it?
### Requirements
 
First, you need to have installed Nextflow (>=20.07) and Singularity.
In order to realize a hac basecalling, you need to have a GPU.

#### Preparation of inputs

#### Run the pipeline

To run run the pipeline, execute:

` nextflow run Nanopore_SARS-CoV-2/main.nf -params-file <params> -profile <profile>`

In <params>, you need to provide inputs, configurations and other options. These are:

| Parameter | Requiered | Default | Description |
| ------ | ------ | --- | ---- |
| run_id | yes | --- | Identificator to the sequencing. |
| fast5 | yes | --- | Path to fast5 directory obtained in the sequencing. |
| fastq | no | --- | Path to fastq directory obtained in the sequencing. Only in case to skip hac basecalling and start pipeline with fastq files. |
| barcodes | yes | --- | Used kit in the sequencing to multiplexed samples. |
| sequencing-summary | yes | --- | Sequencing summary obtained in the sequencing. |
| primers_scheme | yes | nCoV-2019/V3 | Path to primer schemed used in the library preparation. Must be in the actual directory, don't use absolute path. |
| guppy_config | yes | dna_r9.4.1_450bps_hac.cfg | Configuration to used in hac basecalling. More information about available configurations in the guppy manual. |
| guppy_device | no | 'cuda:0,1' | GPU deviced to used in hac basecalling. |
| primers | yes | artic-ncov2019/primer_schemes | Path to directory with all primers scheme. Primer scheme available in artic-ncov2019 repository.  |
| hac_basecalling | yes | True | Boleean value to run hac basecalling. False to start with fastq files, True to star with fast5 files and hac basecalling. |
| len_min_amplicon | yes | 400 | Minimum lengths of the amplicons.  |
| len_max_amplicon | yes | 700 | maximum lengths of the amplicons plus 200.  |
| ourdir | yes | results | Output directory to store the results. |

So, for example, a full execution command with hac basecalling should look like this:

`nextflow run Nanopore_SARS-CoV-2/main.nf --run_id 20210501 --fast5 nanopore_sequencing/fast5/ --barcodes EXP-NBD104 --sequencing-summary nanopore_sequencing/sequencing_summary_FAO03069_83da979a.txt `

And a full execution command without hac basecalling (starting with fastq files) should look like this:

`nextflow run Nanopore_SARS-CoV-2/main.nf --run_id 20210501 --fast5 nanopore_sequencing/fast5/ --fastq nanopore_sequencing/fastq_pass/barcode0[1-6] --barcodes EXP-NBD104 --sequencing-summary nanopore_sequencing/sequencing_summary_FAO03069_83da979a.txt --hac_basecalling False`


Alternatively, you can provide a yaml file containing all the parameters you want to setup (that way you don't have to write everything on the command line). Just download params.example.yml and edit it to your needs (you can delete parameters from the file if you don't want to use them). Then execute the pipeline like this:

` nextflow run Nanopore_SARS-CoV-2/main.nf -params-file params.yml `