# hcov19-nanopore

## ¿Qué incluye?
- Secuencias de consenso para para cada muestra ensamblada.
- Todos los archivos generados por el pipeline ARTIC (bam, fasta, vcf, etc).
- Estadisticas de calidad de las secuncias y ensambles (reportes fastqc, multiqc, archivos de analisis de cobertura y profundidad).
- Reporte que incluye los clades de GISAID, Nextstrain y linaje de Pangolin asociado a cada muestra ensamblada.
- Raw data demultiplexada y filtrada.
- Archivo de GISAID utilizado para la subida de datos en la misma plataforma. Incluye información básica como tipo de virus, fecha de obtención de la muestra, tecnología de secuenciación utilizada, metodo de ensamble y cobterura de profundida obtenida.
- Archivo de resumen qu contiene la información básica de las muestras (en caso de que hayan sido provistas en el archivo samples_data). Estas son Rut, Nombre, Fecha de toma de muestra, Clado de GISAID, Nextstrain, Linaje Pangolin, Porcentaje de cobertura, Media de profundidad del genoma.


## ¿Cómo usarlo?
### Requerimientos

Primero, se necesita tener instalado Nextflow (>=20.07) y Singularity.
En caso de querer realizar el basecalling de alta presición, necesitará contar con una GPU.
Se debe contar con acceso a un directorio que contenga todos los esquemas de primers (disponibles en el repositorio de ARTIC)
#### Preparación de los archivos de entrada

Se debe proveer un archivo donde se describan las muestras, este debe indicarse en el archivo params (sample_data). El archivo debe contener las columnas sample y barcode.
```
sample, barcode
hCoV-19/Chile/MA-CADIUMAG-001/2021, barcode01
hCoV-19/Chile/MA-CADIUMAG-002/2021, barcode02
hCoV-19/Chile/MA-CADIUMAG-003/2021, barcode03
```
#### Ejecución del flujo de trabajo 

Para correr el pipeline ejecutar: 

` nextflow run Nanopore_SARS-CoV-2/main.nf -params-file  <params> -profile <profile>`

En el archivo <params>, se deben proporcionar los archivos de entrada, configuraciones y otras opciones. Estas son:

| Parámetro | Requerido | Por defecto | Descripción |
| ------ | ------ | --- | ---- |
| samples_data | si | --- | Archivo delimitado por comas que debe contener al menos las columnas: sample y barcode, utilizado para definir los nombres de las muestras en procesos posteriores. |
| fast5_directory | si | --- | Ruta al directorio FAST5 obtenido en la secuenciación. |
| fastq_directory | no | --- | Ruta al directorio FASTQ obtenido en la secuenciación en caso de que se haya realizado basecalling, con lo cual se usarán directamente estos datos.  |
| sequencing-summary | no | --- |  Resumen obtenido en el proceso de basecalling. Requerido si se proveen datos en FASTQ. |
| guppy_basecalling_config | no | dna_r9.4.1_450bps_hac.cfg | Configuración a utilizar en el basecalling de alta precisión. Más información sobre las configuraciones disponibles en el manual de Guppy. |
| guppy_barcodes | no | EXP-NBD104 EXP-NBD114  | Kit(s) utilizado(s) en la secuenciación para el multiplexado. |
| guppy_demu_both_ends | no | true | Valor lógico que indica si se realizará la demultiplexación con un solo barcode (false) o dos barcodes (true).   |
| guppy_device | no | auto | Dispositivo GPU a utilizar en el basecalling. |
| artic_primer_schemes_directory | no | --- | Ruta al directorio que contiene todos los esquemas de primers. Esquemas de primers disponibles en el repositorio artic-ncov2019.  |
| artic_primers_scheme | no | nCoV-2019/V3 | Ruta al directorio con el esquema de primers utilizado en la preparación de librería.  |
| artic_normalise | no | 200 |  Establece un valor de cobertura objetivo para reducir los tiempos de ejecución.  |
| gisaid_clades | no | data/gisaid_clades.csv | Ruta al archivo que contiene la descripción de los clados (disponible en repositorio).  |
| gisaid_template | no | data/gisaid_template.csv | Ruta al template utilizado para subir muestras a GISAID (disponible en repositorio).  |
| publish_minimum_coverage | 90 | data/gisaid_template.csv | Valor entre 0 - 100 que indica el porcentaje de bases no identificadas que se permitirán en el consenso para su inclusión en el archivo para subir a GISAID. |  
| output_directory | no | results | Directorio en el cual se almacenaran los resultados. |



Por ejemplo, una ejecución que incluya el basecalling de alta precisión y utilizando los parámtros por defecto, sería algo así:

`nextflow run hcov19-nanopore/main.nf --fast5 input/fast5/  --sample_data input/samples.csv --artic_primer_schemes_directory input/primer_schemes/`

Y una ejecución sin basecalling de alta precisión (comenzando con los archivos fastq demultiplexados), luciria así:
 
`nextflow run hcov19-nanopore/main.nf --fast5 input/fast5/ --sample_data input/samples.csv--fastq input/fastq_pass/barcode0[1-6] --sequencing-summary input/sequencing_summary.txt`

También se puede proporcionar un archivo yaml que contenga todos los parametros a configurar, de esta manera no es necesario escribir todo en la línea de comando. Para esto descargue el archivo de ejemplo params.example.yaml y editelo de acuerdo a sus configuraciones. Luego puede ejecutar el pipeline de la siguiente manera: 

` nextflow run hcov19-nanopore/main.nf -params-file params.yml `



## Resultados 
Una vez terminada la ejecución del pipeline se podrán encontrar los resultados en la carpeta de output definida (results por defecto)
- `artic`: Contiene los resultados asociados a la ejecución del pipeline de ARTIC para el ensamblaje. Dentro de este directorio se pueden encontrar:
	- `bam`:  Contiene los archivos BAM utilizados para la etapa final del ensamblaje para cada una de las muestras.
	- `consensus`: Contiene las secuencias de consenso para cada una de las muestras en formato FASTA.
	- `vcf`: Contiene los archivos VCF (variaciones en la secuencia respecto a la referencia).
	- `pipeline`: Contiene las carpetas completas entregadas en la ejecución del pipeline de ARTIC para cada una de las muestras.
- `qc`:  Contiene todas las estadísticas obtenidas. Dentro de este directorio se pueden encontrar:
	- `alignment_stats/:` Contiene las estadísticas obtenidas con la herramienta samtools de cobertura y profundidad provenientes de los archivos BAM.
	- `fastqc/`: Contiene los archivos de análisis de calidad de las secuencias utilizadas para los consensos de cada muestra secuenciada.
	- `multiqc/`: Contiene un único reporte que incluye todas las métricas de calidad de cada de una de las muestras.
	
- `lineages`: Contiene los archivos en formato csv sobre el genotipado asignado a cada muestra ensamblada, tanto para los clades de GISAID, Nextstrain y linajes de Pangolin
- `raw_data`: Contiene la raw data proveniente de la secuenciación demultiplexada y filtrada por tamaños de lectura válidos de acuerdo a los valores esperados.
- `EpiCoV_BulkUpload.xls`: Archivo en formato Excel definido por GISAID para la subida de datos. Este archivo contiene información básica de las muestras ensambladas como tipo del virus, fecha de obtención de muestra, tecnología de secuenciación utilizada, método de ensamble, cobertura. 
- `samples_summary.xls`: Archivo e formato Excel que contiene la información de cada muestra ensamblada, ésta es:
	- Identificador (sólo si se encontraba en el archivo de entrada)
	- Nombre paciente (sólo si se encontraba en el archivo de entrada)
	- Fecha de toma de muestra (sólo si se encontraba en el archivo de entrada)
	- Clado de GISAID asignado a la muestra
	- Clado de Nextstrain asignado a la muestra
	- Linaje de Pangolin asignado a la muestra
	- Porcentaje de cobertura de la referencia
	- Media de la profundidad del genoma
	

## What is included?
- Consensus sequences for each assembled sample.
- All files generated by the ARTIC pipeline  (bam, fasta, vcf, etc)
- Quality statistics of the sections and assemblies  (fastqc and multiqc reports, and,  coverage, depth analysis files)
- Report that includes the clades of GISAID, Nextstrain and lineage of Pangolin associated with each assembled sample. 
- Raw data demultiplexed y filtered.
- GISAID file used for uploading data on the same platform. It includes basic information such as type of virus, date of obtaining the sample, sequencing technology used, assembly method and coverage of depth obtained. .
- Summary file containing the basic information of the samples (in case they have been provided in the samples_data file). These are Rut, Name, Date of sample collection, GISAID Clado, Nextstrain, Pangolin lineage, Percentage of coverage, Average depth of the genome. 

## How to use it?
### Requirements
 
First, you need to have installed Nextflow (>=20.07) and Singularity.
In order to realize a hac basecalling, you need to have a GPU.
You must have access to a directory that contains all the primer schemes (available in the ARTIC repository) 
#### Preparation of inputs
A file must be provided where the samples are described, this must be indicated in the params file (sample_data). The file must contain at least the columns sample and barcode.
```
sample, barcode
hCoV-19/Chile/MA-CADIUMAG-001/2021, barcode01
hCoV-19/Chile/MA-CADIUMAG-001/2021, barcode02
hCoV-19/Chile/MA-CADIUMAG-001/2021, barcode03
```
#### Run the pipeline

To run run the pipeline, execute:

` nextflow run Nanopore_SARS-CoV-2/main.nf -params-file <params> -profile <profile>`

In <params>, you need to provide inputs, configurations and other options. These are:

| Parameter | Requiered | Default | Description |
| ------ | ------ | --- | ---- |
| samples_data | yes | --- | Comma delimited file that must contain at least the columns: sample and barcode, used to define the names of the samples in subsequent processes.  |
| fast5_directory | yes | --- | Path to fast5 directory obtained in the sequencing. |
| fastq_directory | no | --- | Path to fastq directory obtained in the sequencing if the basecalling was carried out. If the directory is provided, the high precision basecalling will not be performed.  |
| sequencing-summary | no | --- | Sequencing summary obtained in the basecalling.Required if data is provided in FASTQ.  |
| guppy_basecalling_config | no | dna_r9.4.1_450bps_hac.cfg | Configuration to used in hac basecalling. More information about available configurations in the guppy manual. |
| guppy_barcodes | no | EXP-NBD104 EXP-NBD114   | Used kit in the sequencing to multiplexed samples. |
| guppy_demu_both_ends | no | true | Logical value that indicates whether demultiplexing will be performed with a single barcode (false) or two barcodes (true).    |
| guppy_device | no | auto | GPU deviced to used in hac basecalling. |
| artic_primer_schemes_directory | no |  | Path to directory with all primers scheme. Primer scheme available in artic-ncov2019 repository.  |
| artic_primers_scheme | no | nCoV-2019/V3 | Path to primer schemed used in the library preparation. Must be in the actual directory, don't use absolute path. |
| artic_normalise | no | 200 |  Set a target coverage value to reduce execution times.   |
| gisaid_clades | no | data/gisaid_clades.csv | Path to the file that contains the description of the clades (available in the repository).  |
| gisaid_template | no | data/gisaid_template.csv | Path to the template used to upload samples to GISAID (available in repository).  |
| publish_minimum_coverage | 90 | data/gisaid_template.csv | Value between 0 - 100 that indicates the percentage of unidentified bases that will be allowed in the consensus for inclusion in the file to upload to GISAID.  |  
| output_directory | no | results | Output directory to store the results. |


So, for example, a full execution command with hac basecalling should look like this:

`nextflow run hcov19-nanopore/main.nf --fast5 input/fast5/  --sample_data input/samples.csv --artic_primer_schemes_directory input/primer_schemes/`

And a full execution command without hac basecalling (starting with fastq files) should look like this:

`nextflow run hcov19-nanopore/main.nf --fast5 input/fast5/ --sample_data input/samples.csv--fastq input/fastq_pass/barcode0[1-6] --sequencing-summary input/sequencing_summary.txt`


Alternatively, you can provide a yaml file containing all the parameters you want to setup (that way you don't have to write everything on the command line). Just download params.example.yml and edit it to your needs (you can delete parameters from the file if you don't want to use them). Then execute the pipeline like this:

` nextflow run hcov19-nanopore/main.nf -params-file params.yml `


## Results 
Once the pipeline finished running you will find a set of files. These are:
- `artic`: Contains the results associated with the execution of the ARTIC pipeline for the assembly. Within this directory you can find: 
	- `bam`:  Contains the BAM files used for the final stage of assembly for each of the samples. 
	- `consensus`: 	Contains the consensus sequences for each of the samples in FASTA format. 
	- `vcf`: Contains the VCF files (variations in sequence from reference). 
	- `pipeline`: Contains the complete folders delivered in the execution of the ARTIC pipeline for each of the samples. 
- `qc`:  Contiene todas las estadísticas obtenidas. Dentro de este directorio se pueden encontrar:
	- `alignment_stats/:` Contains the statistics obtained with the coverage and depth samtools tool from the BAM files. 
	- `fastqc/`: Contains the quality analysis files of the sequences used for the consensuses of each sequenced sample. 
	- `multiqc/`: Contains a single report that includes all the quality metrics for each of the samples. 
	
- `lineages`: Contains the files in csv format on the genotyping assigned to each assembled sample, for both the GISAID, Nextstrain and Pangolin lineages. 
- `raw_data`: Contains the raw data from the demultiplexed sequencing and filtered by valid read sizes according to the expected values. 
- `EpiCoV_BulkUpload.xls`: File in Excel format defined by GISAID for uploading data. This file contains basic information on the assembled samples such as virus type, sample collection date, sequencing technology used, assembly method, coverage. 
- `samples_summary.xls`: File in Excel format that contains the information of each assembled sample, this is: 
	- Identifier (only if it was in the input file )
	- Patient name  (only if it was in the input file a)
	- Sample collection date  (only if it was in the input file )
	- GISAID clade assigned to the sample 
	- Nextstrain clade assigned to sample 
	- Pangolin lineage assigned to the sample 
	- Reference coverage percentage 
	- Genome depth mean
