# Pipeline para Análisis de nCoV-2019 con Oxford Nanopore

[![DOI](https://zenodo.org/badge/367064011.svg)](https://zenodo.org/badge/latestdoi/367064011)

Pipeline en Nextflow para el ensamblaje y posterior análisis de genomas de SARS-CoV-2.
Utiliza el protocolo bioinformático de ARTIC, por lo que es necesario haber utilizado el protocolo de ARTIC para la secuenciación ([más información](https://artic.network/ncov-2019)).

Creado para sistematizar el análisis bioinformático asociado a la vigilancia genómica realizada en el Centro Asistencial Docente y de Investigación de la Universidad de Magallanes. Más información sobre nosotros en [catg.cl](https://catg.cl).

Actualmente incluye:

- Basecalling y demultiplexado usando Guppy (opcional, se requiere Guppy instalado)
- Ensamblaje y obtención de variantes usando [ARTIC](https://github.com/artic-network/fieldbioinformatics)
- Anotación de variantes usando [SnpEff](https://pcingola.github.io/SnpEff/) (gen, cambio de aminoácido)
- Identificación de linajes de GISAID, Pangolin y Nexstrain
- Recolección de métricas de calidad y estadísticas (qc, coverage, etc)
- Generación de resumen final en formato Excel incluyendo la información más relevante obtenida en la ejecución del pipeline
- Preparación de planilla y secuencias para ser subidas a GISAID

## Inicio Rápido

1. Instalar Nextflow, Singularity (o Docker) y Guppy (si se desea hacer basecalling)
2. Crear un archivo CSV con las columnas `barcode` y `sample`, que tenga los barcodes de Nanopore (barcode01, barcode02, etc) y un nombre de muestra para cada uno de ellos. Si se planea la carga en GISAID, es recomendable tener una tercera columna: `gisaid_covv_virus_name`, que tenga los nombres que se utilizarán para subir las muestras a GISAID (debe tener el formato `hCoV-19/<Pais>/<ID>/<Año>`).
3. Ejecutar el pipeline. Ejemplos:

   Con basecalling, singularity

   ```
   nextflow run catg-umag/ncov2019-ont-nf -r v1.0 -profile singularity \
       --sample_data <csv_file> --fast5_directory <fast5_dir>
   ```

   Sin basecalling, docker

   ```
   nextflow run catg-umag/ncov2019-ont-nf -r v1.0 -profile docker \
       --sample_data <csv_file> --fast5_directory <fast5_dir> \
       --fastq_directory <fastq_dir> --sequencing_summary <txt_file>
   ```

4. Esperar a que finalice la ejecución. En `results/` se encontrarán los archivos entregados por el pipeline.

## Requisitos

- [Nextflow](https://www.nextflow.io/) (>= 20.07)
- [Singularity](https://sylabs.io/guides/3.7/admin-guide/) o [Docker](https://www.docker.com/get-started)
- Guppy (sólo si se desea hacer basecalling, se descarga desde la página de la [comunidad de Oxford Nanopore](https://community.nanoporetech.com/downloads))

Debido a Nextflow, los sistemas operativos soportados son GNU/Linux y MacOS. En caso de querer de todas maneras utilizarlo en Windows, es necesario usar [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10).

## Utilización del Pipeline

### Entradas

Dependiendo del punto de partida existen dos casos:

1. Ya se hizo basecalling y demultiplexing ya sea en tiempo real durante la secuenciación u otra instancia, y no se pretende hacer nuevamente. En este se requerirá:
   - el sub-directorio `pass` del directorio `fast5`
   - el directorio con el resultado del basecalling + demultiplexing; en este directorio deberían encontrarse los directorios `barcodeXX` correspondientes a cada muestra
   - el archivo `sequencing_summary.txt` que puede estar dentro del directorio de los datos crudos o con basecalling
   - un archivo CSV con el detalle de las muestras
2. Se desea hacer basecalling como parte del pipeline, ya sea porque no se hizo previamente o quiere realizarse basecalling de alta precisión. Hay que tener en consideración que el basecalling de alta precisión en considerablemente más intensivo computacionalmente, por lo que se hace casi indispensable hacerlo utilizando una GPU soportada. En este caso, las entradas son:
   - el directorio `fast5` en su totalidad
   - el archivo CSV con el detalle de las muestras

### Cómo preparar el archivo de detalle de muestras

Este archivo (en formato CSV) debe contar con al menos dos columnas: `barcode` y `sample`, en el cual se asocie cada barcode a una muestra determinada.
El campo sample será utilizado para nombrar todos los archivos asociados a cada muestra, por lo que es ideal que sea algo sencillo, ojalá sólo letras y números. Además de las dos columnas requeridas se pueden incluir columnas adicionales, las cuales serán "traspasadas" al resumen final.

Otro conjunto de campos que puede incluir este archivo son relevantes si se desea subir los datos a GISAID. Estos corresponden a cualquiera de las columnas de la plantilla de carga de datos de GISAID (disponible en `data/`), que deben ser indicados con el prefijo `gisaid_` seguido del nombre en la primera fila de la plantilla. Como recomendación, sería conveniente incluir al menos `covv_virs_name`, ya que este campo será utilizado para renombrar las secuencias en el archivo FASTA de forma de que correspondan a los nombres indicados en esta columna, proceso que deberá realizarse de forma manual en caso de no proveer estos valores. Revisar la plantilla para ver cuáles son los campos disponibles y cómo deben ser llenados. Algunas columnas serán llenadas aunque no sean especificadas, éstas son: `fn`, `covv_type`, `covv_passage`, `covv_host`, `covv_seq_technology`, `covv_assembly_method` y `covv_coverage`. Será necesario rellenar manualmente los campos obligatorios que falten para que la plataforma acepte la carga.

Por ejemplo, para las siguientes 5 muestras se incluyó el campo adicional `city` (el cual figurará en el resumen final), además de `covv_virus_name`, `covv_collection_date` y `covv_location` para ser utilizados en la plantilla de carga.

```
sample,barcode,city,gisaid_covv_virus_name,gisaid_covv_collection_date,gisaid_covv_location
2101011,barcode01,Punta Arenas,hCoV-19/Chile/CADIUMAG-51/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2101018,barcode02,Punta Arenas,hCoV-19/Chile/CADIUMAG-52/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2101053,barcode03,Punta Arenas,hCoV-19/Chile/CADIUMAG-53/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2101025,barcode04,Punta Arenas,hCoV-19/Chile/CADIUMAG-54/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
2105001,barcode05,Punta Arenas,hCoV-19/Chile/CADIUMAG-55/2021,2021-05-20,South America / Chile / Magallanes / Punta Arenas
```

### Parámetros

Además de la posibilidad de especificar los inputs, es posible controlar algunos aspectos del pipeline a través de parámetros. El listado de parámetros disponibles es el siguiente:

| Parámetro                      | Requerido | Valor Por defecto            | Descripción                                                                                                                                                    |
| ------------------------------ | --------- | ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| sample_data                    | si        | ---                          | Archivo delimitado por comas con el detalle de las muestras.                                                                                                   |
| fast5_directory                | si        | ---                          | Ruta al directorio FAST5 obtenido en la secuenciación. En caso no realizar basecalling, indicar el sub-directorio `pass`.                                      |
| fastq_directory                | no        | ---                          | Ruta al directorio FASTQ obtenido en la secuenciación en caso de que se haya realizado basecalling y no se desee realizar basecalling.                         |
| sequencing_summary             | no        | ---                          | Resumen obtenido en el proceso de basecalling. Requerido si se proveen datos en FASTQ.                                                                         |
| run_id                         | no        | ""                           | ID opcional para la secuenciación que se utilizará como sufijo en los nombres de archivos de los resúmenes y recopilaciones de datos.                          |
| guppy_basecalling_config       | no        | dna_r9.4.1_450bps_hac.cfg    | Configuración a utilizar para realizar el basecalling con guppy. Más información sobre las configuraciones disponibles en el manual de Guppy.                  |
| guppy_barcodes                 | no        | "EXP-NBD104 EXP-NBD114"      | Kit(s) utilizado(s) en la secuenciación para la multiplexación.                                                                                                |
| guppy_demux_both_ends          | no        | true                         | Valor que indica si se debe exigir la presencia de ambos barcodes (5' y 3') al momento de hacer demultiplexing.                                                |
| guppy_cpus                     | no        | 16                           | Cantidad de CPUs a utilizar en los procesos asociados a Guppy.                                                                                                 |
| guppy_basecalling_extra_config | no        | "--device auto"              | Opciones extras para Guppy al hacer basecalling (por ejemplo, parámetros asociados a la configuración de la GPU).                                              |
| artic_primers_scheme           | no        | nCoV-2019/V3                 | Esquema de primers de ARTIC utilizado al construir la librería.                                                                                                |
| artic_normalise                | no        | 500                          | Valor para establecer un valor de cobertura objetivo en el pipeline de ARTIC.                                                                                  |
| gisaid_clades                  | no        | data/gisaid_clades.csv       | Ruta al archivo que contiene la descripción de los clados de GISAID (disponible en repositorio).                                                               |
| gisaid_template                | no        | data/20210222_EpiCoV....xlsx | Ruta al template utilizado para subir muestras a GISAID (disponible en repositorio).                                                                           |
| gisaid_submission_enabled      | no        | true                         | Habilitar (o no) la generación de los archivos preparados para cargar datos a GISAID.                                                                          |
| publish_minimum_completion     | no        | 95                           | Valor entre 0 - 100 que indica el porcentaje de bases cubiertas respecto al genoma que se requerirá para la inclusión de estas muestras en la carga de GISAID. |
| output_directory               | no        | results                      | Directorio en el cual se almacenaran los resultados.                                                                                                           |

El listado y lo valores por defecto se encuentran en el archivo `params.default.yml`.
Los parámetros pueden ser indicados a través de la línea de comandos (por ejemplo `--run_id 20210501A`), pero también pueden pasarse a través de un archivo en formato YAML, se puede utilizar como plantilla el archivo `params.default.yml` (no mover o editar el archivo en sí, ya que es necesario).

### Ejecución

El pipeline puede descargarse directamente y e indicar a Nextflow la ruta donde este se encuentra, pero también puede indicarse directamente `catg-umag/ncov2019-ont-n` y una versión a través del argumento `-r`, lo que descargará automáticamente el pipeline. Por ejemplo:

```
nextflow run catg-umag/ncov2019-ont-nf -r v1.0 -profile singularity ...
```

Es imperativo utilizar el perfil `singularity` o `docker` para que utilice contenedores para acceder a las herramientas requeridas.

Ejemplo con basecalling utilizando parámetros a través de línea de comandoss:

```
nextflow run ncov2019-ont-nf/ -profile singularity --fast5_directory input/fast5/ \
    --sample_data input/samples.csv --run_id R210505
```

Ejemplo utilizando parámetros mediant archivo YAML:

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

El pipeline por defecto se ejecutará de manera local, pero en caso de querer ejecutarlo en un clúster de cómputo que utilice Slurm, puede especificarse mediante el perfil `slurm`. Y no olvidar también el perfil de singularity (este debe estar instalado en el clúster): `-profile slurm,singularity`.

## Resultados

Dentro del directorio de resultados (`results/` por defecto), se encontrará lo siguiente:

- `raw_data`: archivos FASTQ para cada muestra, generado después de los procesos de basecalling, demultiplexing y filtro de largo
- `artic/`: resultados de la ejecución del pipeline de ARTIC para cada una de las muestras
- `qc/`: métricas de calidad y cobertura del genoma para cada una de las muestras
- `vcf/`: archivos VCF (variaciones en el genoma) para cada una de las muestras, anotadas con SnpEff
- `lineages/`: resultados sin procesar de las identificaciones de linajes
- `summary/`: directorio con resúmenes de la información más relevante recolectada en el pipeline:
  - `all_consensus.fasta`: archivo FASTA con los consensos para todas las muestras (incluso las incompletas)
  - `sample_summary.csv`: archivo CSV con la información incluída en el archivo de entrada, linajes y métricas de cobertura
  - `variants_list.csv`: archivo CSV con las variantes para todas las muestras con su información esencial
  - `summary.xlsx`: archivo Excel con la información de las muestras en la primera hoja, y un resúmen de presencia de las diferentes variantes en las muestras en la segunda (las cuales son reordenadas dependiendo de su similitud considerando las variantes que poseen).
- `gisaid_submission`: si fue habilitada la opción de preparar la carga a GISAID, en este directory se encontrará un archivo FASTA filtrado por el criterio de cobertura y con las secuencias renombrados si se indicaron los nombres, además de la planilla Excel llenado con la información disponible (podría aún requerir llenar campos manualmente)

## Cómo Citar

Si este trabajo te fue de utilidad, puedes citarlo a través del siguiente artículo:

> González-Puelma, J.; Aldridge, J.; Montes de Oca, M.; Pinto, M.; Uribe-Paredes, R.; Fernández-Goycoolea, J.; Alvarez-Saravia, D.; Álvarez, H.; Encina, G.; Weitzel, T.; Muñoz, R.; Olivera-Nappa, Á.; Pantano, S.; Navarrete, M.A. Mutation in a SARS-CoV-2 Haplotype from Sub-Antarctic Chile Reveals New Insights into the Spike’s Dynamics. Viruses 2021, 13, 883. https://doi.org/10.3390/v13050883

O también puedes usar la [referencia de Zenodo](https://zenodo.org/badge/latestdoi/367064011).

## Software Utilizado

- ARTIC y todas sus dependencias
- Python y sus librerías: biopython, openpyxl, pandas, scipy
- Herramientas bioinformáticas: FastQC, MultiQC, samtools, SnpEff, Nextclade, Pangolin
- Otros: Nextflow, unoconv
