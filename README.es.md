# Pipeline para Análisis de nCoV-2019 con Oxford Nanopore

Pipeline en Nextflow para el ensamblaje y posterior análisis de genomas de SARS-CoV-2.
Utiliza el protocolo bioinformático de ARTIC, por lo que es necesario haber utilizado el protocolo de ARTIC para la secuenciación ([más información](https://artic.network/ncov-2019)).

Desarrollado para sistematizar el análisis bioinformático asociado a la vigilancia genómica realizada en el Centro Asistencial Docente y de Investigación de la Universidad de Magallanes. Más información sobre nosotros en [catg.cl](https://catg.cl).

Actualmente incluye:

- Ensamblaje y obtención de variantes usando [ARTIC](https://github.com/artic-network/fieldbioinformatics)
- Anotación de variantes usando [SnpEff](https://pcingola.github.io/SnpEff/) (gen, cambio de aminoácido)
- Identificación de linajes de Pangolin y Nexstrain (incluyendo nombre común de variante)
- Recolección de métricas de calidad y estadísticas (qc, coverage, etc)
- Generación de resumen final en formato Excel incluyendo la información más relevante obtenida en la ejecución del pipeline
- Preparación de planilla y secuencias para ser subidas a GISAID

## Inicio Rápido

1. Instalar Nextflow y Singularity (o Docker)
2. Crear un archivo CSV con las columnas `sample` y `fastq_file`, que asocie los nombres de las muestras con los datos en FASTQ. Si se planea la carga en GISAID, es recomendable tener una tercera columna: `gisaid_covv_virus_name`, que tenga los nombres que se utilizarán para subir las muestras a GISAID (debe tener el formato `hCoV-19/<Pais>/<ID>/<Año>`).
3. Ejecutar el pipeline. Ejemplos:

   Usando Apptainer:

   ```bash
   nextflow run catg-umag/ncov2019-ont-nf -r v2.0.0 -profile apptainer \
       --sample_data <csv_file>
   ```

    Usando Nanopolish en lugar de Medaka, con Docker:

   ```bash
   nextflow run catg-umag/ncov2019-ont-nf -r v2.0.0 -profile docker \
       --sample_data <csv_file> --artic_use_nanopolish \
       --fast5_directory <fast5_dir> --sequencing_summary <txt_file>
   ```

4. Esperar a que finalice la ejecución. En `results/` se encontrarán los archivos entregados por el pipeline.

## Requisitos

- [Nextflow](https://www.nextflow.io/) (>= 22.04.0)
- [Apptainer](https://apptainer.org/docs/admin/main/installation.html) o [Docker](https://www.docker.com/get-started)

Debido a Nextflow, los sistemas operativos soportados son GNU/Linux y MacOS. En caso de querer de todas maneras utilizarlo en Windows, es necesario usar [WSL](https://learn.microsoft.com/en-us/windows/wsl/install).

## Utilización del Pipeline

### Entradas

Este pipeline está diseñado para trabajar con datos después de haber realizado basecalling y demultiplexación, o sea los archivos FASTQ. De no tener ese proceso realizado, o querer realizar basecalling con mayor calidad, se puede realizar con [este pipeline](https://github.com/catg-umag/ont-basecalling-demultiplexing).

Entonces, las entradas requeridas son:
- Un archivo CSV con el detalle de las muestras (ver más abajo)
- Un archivo FASTQ único por muestra (comprimido o no). Guppy por defecto genera una carpeta con múltiples archivos FASTQ, por lo que si se desea utilizar estos datos, es necesario concatenarlos antes.
- (Sólo si se desea utilizar Nanopolish en lugar de Medaka). El subdirectorio `pass` del directorio con los archivos FASTQ generados el basecalling, y también el archivo `sequencing_summary.txt`

### Cómo preparar el archivo de detalle de muestras

Este archivo (en formato CSV) debe contar con al menos dos columnas: `sample` y `fastq_file`, en el cual se asocie cada muestra a un archivo FASTQ.
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

| Parámetro                  | Requerido | Valor Por defecto            | Descripción                                                                                                                                                    |
| -------------------------- | --------- | ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| sample_data                | sí        | ---                          | Archivo delimitado por comas con el detalle de las muestras.                                                                                                   |
| fast5_directory            | no        | ---                          | Ruta al subdirectorio `pass` de FAST5 obtenido en la secuenciación, sólo en caso de ocupar Nanopolish.                                                         |
| sequencing_summary         | no        | ---                          | Resumen obtenido en el proceso de basecalling. Requerido si se proveen datos en FASTQ.                                                                         |
| run_id                     | no        | ""                           | ID opcional para la secuenciación que se utilizará como sufijo en los nombres de archivos de los resúmenes y recopilaciones de datos.                          |
| artic_primer_scheme        | no        | SARS-CoV-2/V4.1              | Esquema de primers de ARTIC utilizado al construir la librería.                                                                                                |
| artic_medaka_model         | no        | r103_hac_g507                | Modelo de Medaka (debe coincidir con el modelo utilizado par el basecalling).                                                                                  |
| artic_use_nanopolish       | no        | false                        | Usar Nanopolish en lugar de Medaka para el pulido.                                                                                                             |
| artic_normalise            | no        | 500                          | Valor para establecer un valor de cobertura objetivo en el pipeline de ARTIC.                                                                                  |
| gisaid_template            | no        | data/20210222_EpiCoV....xlsx | Ruta al template utilizado para subir muestras a GISAID (disponible en repositorio).                                                                           |
| gisaid_submission_enabled  | no        | true                         | Habilitar (o no) la generación de los archivos preparados para cargar datos a GISAID.                                                                          |
| publish_minimum_completion | no        | 95                           | Valor entre 0 - 100 que indica el porcentaje de bases cubiertas respecto al genoma que se requerirá para la inclusión de estas muestras en la carga de GISAID. |
| multiqc_config             | no        | conf/multiqc_config.yaml     | Archivo de configuración de MultiQC (ya se incluye uno).                                                                                                       |
| output_directory           | no        | results                      | Directorio en el cual se almacenaran los resultados.                                                                                                           |

El listado y lo valores por defecto se encuentran en el archivo `params.default.yml`.
Los parámetros pueden ser indicados a través de la línea de comandos (por ejemplo `--run_id 20210501A`), pero también pueden pasarse a través de un archivo en formato YAML, se puede utilizar como plantilla el archivo `params.default.yml` (no mover o editar el archivo en sí, ya que es necesario).

### Ejecución

El pipeline puede descargarse directamente y e indicar a Nextflow la ruta donde este se encuentra, pero también puede indicarse directamente `catg-umag/ncov2019-ont-n` y una versión a través del argumento `-r`, lo que descargará automáticamente el pipeline. Por ejemplo:

```bash
nextflow run catg-umag/ncov2019-ont-nf -r v2.0.0 -profile apptainer ...
```

Es imperativo utilizar el perfil `singularity` o `docker` para que utilice contenedores para acceder a las herramientas requeridas. De lo contrario, no funcionará.

Ejemplo con basecalling utilizando parámetros a través de línea de comandos:

```bash
nextflow run ncov2019-ont-nf/ -profile apptainer --sample_data input/samples.csv --run_id R210505
```

Ejemplo utilizando parámetros mediant archivo YAML:

`params.yml`:

```yaml
sample_data: input/samples.csv
run_id: R210505
artic_medaka_model: r941_min_hac_g507
artic_normalise: 1000
gisaid_submission_enabled: false
```

```bash
nextflow run ncov2019-ont-nf/ -profile apptainer -params-file params.yml
```

El pipeline por defecto se ejecutará de manera local, pero en caso de querer ejecutarlo en un clúster de cómputo que utilice Slurm, puede especificarse mediante el perfil `slurm`. Y no olvidar también el perfil de singularity (este debe estar instalado en el clúster): `-profile slurm,singularity`.

## Resultados

Dentro del directorio de resultados (`results/` por defecto), se encontrará lo siguiente:

- `fastq_data/`: archivos FASTQ para cada muestra, después de filtros de largo
- `artic/`: resultados de la ejecución del pipeline de ARTIC para cada una de las muestras
- `qc/`: métricas de calidad y cobertura del genoma para cada una de las muestras
- `vcf/`: archivos VCF (variaciones en el genoma) para cada una de las muestras, anotadas con SnpEff
- `lineages/`: resultados sin procesar de las identificaciones de linajes
- `summary/`: directorio con resúmenes de la información más relevante recolectada en el pipeline:
  - `all_consensus.fasta`: archivo FASTA con los consensos para todas las muestras (incluso las incompletas)
  - `sample_summary.csv`: archivo CSV con la información incluída en el archivo de entrada, linajes y métricas de cobertura
  - `variants_list.csv`: archivo CSV con las variantes para todas las muestras con su información esencial
  - `summary.xlsx`: archivo Excel con la información de las muestras en la primera hoja, y un resúmen de presencia de las diferentes variantes en las muestras en la segunda (las cuales son reordenadas dependiendo de su similitud considerando las variantes que poseen).
- `gisaid_submission/`: si fue habilitada la opción de preparar la carga a GISAID, en este directory se encontrará un archivo FASTA filtrado por el criterio de cobertura y con las secuencias renombrados si se indicaron los nombres, además de la planilla Excel llenado con la información disponible (podría aún requerir llenar campos manualmente)

## Cómo Citar

Si este trabajo te fue de utilidad, puedes citarlo a través del siguiente artículo:

> González-Puelma, J.; Aldridge, J.; Montes de Oca, M.; Pinto, M.; Uribe-Paredes, R.; Fernández-Goycoolea, J.; Alvarez-Saravia, D.; Álvarez, H.; Encina, G.; Weitzel, T.; Muñoz, R.; Olivera-Nappa, Á.; Pantano, S.; Navarrete, M.A. Mutation in a SARS-CoV-2 Haplotype from Sub-Antarctic Chile Reveals New Insights into the Spike’s Dynamics. Viruses 2021, 13, 883. https://doi.org/10.3390/v13050883

## Software Utilizado

- ARTIC y todas sus dependencias
- Python y sus librerías: biopython, openpyxl, pandas, scipy
- Herramientas bioinformáticas: FastQC, MultiQC, samtools, SnpEff, Nextclade, Pangolin
- Otros: Nextflow, libreoffice
