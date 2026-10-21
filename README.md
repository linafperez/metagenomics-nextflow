# Pipeline de metagenómica shotgun

Pipeline Nextflow DSL2 para reads paired-end locales o una selección explícita de BioSamples de SRA. Reconstruye MAGs con MEGAHIT y metaSPAdes, refina ambos catálogos y produce taxonomía, filogenómica, anotación funcional, abundancia y MultiQC.

## Entrada y ejecución

Cada ejecución utiliza **un solo workflow de Nextflow**. El TSV completo se lee y valida al inicio con `splitCsv` y Groovy. Todas sus muestras entran al dataflow; Nextflow y el executor deciden qué tareas pueden ejecutarse simultáneamente.

```text
TSV de BioSamples + BioProject
       + RunInfo NCBI (una tarea de metadata, o CSV existente)
       -> selección/validación Groovy de TODOS los BioSamples
       -> channel de runs -> SRATOOLS_ACQUIRE por run
       -> reunión de los runs de CADA BioSample
       -> tuple(meta, [R1, R2]) por muestra
       -> FastQC raw -> fastp -> FastQC clean -> Bowtie2 humano
       -> reads no humanas por muestra -> METAGENOMICS_GLOBAL
```

Una muestra puede comenzar FastQC mientras otras siguen descargándose. La reunión de runs utiliza `groupKey` con el número esperado de runs de esa muestra; no espera a las demás. Bowtie2 recibe un par por muestra. Sólo se combinan runs pertenecientes al mismo BioSample. Las coasambleas posteriores conservan su comportamiento científico original y sí necesitan reunir sus muestras.

### Manifest SRA

Se necesitan `--sra-project PRJNA...` y `--sra-samples FILE.tsv` en el launcher, o `--sraProject` y `--sraSamples` al ejecutar Nextflow directamente.

```tsv
biosample_accession	group
SAMN29531760	PD
SAMN29531788	CT
```

- `biosample_accession`: obligatorio, único; admite accesiones `SAMN`, `SAMEA` y `SAMD`.
- `group`: opcional; se detecta automáticamente si existe. `--group-column NAME` permite seleccionar otra columna. Las etiquetas se conservan sin interpretarlas: `CT` es la etiqueta de controles del manifest PRJNA782492; no se convierte en otro nombre.
- Las columnas adicionales se conservan en `meta.metadata` y en `pipeline_info/sample_metadata.tsv`.
- Si se proporciona una columna de grupo, sus valores deben estar completos. Se rechazan identidades duplicadas, columnas duplicadas, campos de más o de menos y metadata contradictoria.

Cada `meta` contiene `id`, `biosample_accession`, `group`, `project_accession`, `run_accessions`, `identity_source`, `single_end` y `metadata`. La posición de una muestra en el TSV no forma parte del hash de adquisición; `sample_order` se escribe en la tabla de metadata del conjunto.

La resolución de RunInfo conserva los criterios existentes: pertenencia al BioProject, acceso público con ruta de descarga, biblioteca `PAIRED`, estrategia `WGS`, origen `METAGENOMIC`, plataforma permitida y contadores consistentes. Por defecto se permiten `ILLUMINA,BGISEQ`. Se conservan todos los runs elegibles de cada BioSample y se ordenan por accesión. Un BioSample sin runs elegibles o un run duplicado causa error antes de descargar reads. Los runs excluidos se registran en `pipeline_info/sra_exclusions.tsv`.

Opcionalmente `--sra-run-info /ruta/RunInfo.csv` (`--sraRunInfo` directamente) utiliza un CSV NCBI existente sin consultar metadata por internet. Para reproducibilidad, conservar ese CSV congelado junto con su fecha, procedencia y checksum. Sigue aplicando la misma validación de elegibilidad. La descarga de reads real requiere después conectividad SRA.

### FASTQ locales

`--input samples.csv` acepta:

```csv
sample,fastq_1,fastq_2,group
sample_A,reads/A_R1.fastq.gz,reads/A_R2.fastq.gz,PD
```

Las rutas relativas se resuelven respecto al CSV original. Se validan archivos, extensiones, identidades únicas y mates distintos con Groovy/Nextflow. La entrada local y SRA son mutuamente excluyentes y comparten el mismo preprocesamiento y análisis posterior.

## Ocho muestras y conjunto completo

Los manifests incluidos son:

| Archivo | BioSamples | Grupos |
| --- | ---: | --- |
| `assets/raw/PRJNA782492_subset8.tsv` | 8 | 4 PD + 4 CT |
| `assets/raw/PRJNA782492_metagenomic_samples.tsv` | 95 | Etiquetas originales PD/CT |

El conjunto completo contiene **95**, confirmado durante el refactor; no se inventó una muestra número 96. El pipeline admite cualquier cantidad de filas. La única selección del conjunto es el manifest, sin ramas de código para ocho muestras.

Con Nextflow >=26.04.6, un runtime disponible y las bases de datos ya preparadas, establecer rutas reales:

```bash
DB_CONFIG=/ruta/real/metagenomics_databases.config
RUN_ROOT=/ruta/persistente/PRJNA782492
```

Validación inicial, ocho muestras:

```bash
bash metagenomics_pipeline.sh --hpc --apptainer --run \
  --sra-project PRJNA782492 \
  --sra-samples assets/raw/PRJNA782492_subset8.tsv \
  --group-column group --database-config "$DB_CONFIG" \
  --outdir "$RUN_ROOT/results" --work-dir "$RUN_ROOT/work"
```

Conjunto completo, mismo código y opciones:

```bash
bash metagenomics_pipeline.sh --hpc --apptainer --run \
  --sra-project PRJNA782492 \
  --sra-samples assets/raw/PRJNA782492_metagenomic_samples.tsv \
  --group-column group --database-config "$DB_CONFIG" \
  --outdir "$RUN_ROOT/results" --work-dir "$RUN_ROOT/work"
```

Para ampliar una ejecución previa, añadir `--resume` y conservar directorio de lanzamiento, `work/` y caché. Para comparar resultados independientes, utilizar otro `RUN_ROOT`. Una ampliación recalcula las coasambleas y análisis dependientes del conjunto; las tareas por muestra cuyos inputs y metadata permanecen iguales pueden reutilizarse.

También se puede ejecutar sin launcher:

```bash
nextflow -c "$DB_CONFIG" run main.nf -profile hpc,apptainer \
  --sraProject PRJNA782492 --sraSamples assets/raw/PRJNA782492_subset8.tsv \
  --groupColumn group --outdir "$RUN_ROOT/results" -work-dir "$RUN_ROOT/work"
```

## Preprocesamiento y decisiones científicas

Se mantienen los módulos y argumentos científicos existentes:

1. SRA Toolkit 3.4.1: `prefetch`, `vdb-validate`, `fasterq-dump --split-files --skip-technical`; validación completa de estructura FASTQ, IDs, mates y longitudes mediante AWK; compresión y comprobación CRC con pigz 2.8. Los singleton reads se descartan, como en la entrada paired-end anterior.
2. FastQC 0.12.1 sobre reads raw. fastp 1.0.1 conserva Q30, longitud mínima 30 nt, detección de adaptadores paired-end, recorte poly-G/poly-X y ventana derecha de cuatro bases con calidad media Q30. FastQC vuelve a evaluar las reads limpias.
3. Bowtie2 2.5.4 conserva `--very-sensitive-local --phred33` y genera las reads no humanas de **cada muestra**, usando el índice humano externo.
4. MEGAHIT y metaSPAdes mantienen sus coasambleas y MetaQUAST. CoverM conserva cobertura y BAMs; COMEBin, MetaBAT2, SemiBin2 y Vamb alimentan DAS Tool.
5. CheckM2 y GUNC evalúan los bins. La selección HQ conserva completitud >90% y contaminación <5%. dRep conserva 99% ANI, agrupación de especies al 95% y cobertura mínima 30%; se integra y evalúa el catálogo final.
6. Se conservan GTDB-Tk, PhyloPhlAn con IQ-TREE, GeneMarkS-2, eggNOG-mapper, InterProScan, integración funcional y abundancia CoverM. La tabla larga conserva grupos y excluye la fila especial `unmapped` sin renormalizar MAGs.
7. MultiQC 1.35 conserva `multiqc_output/multiqc_report_data` y publica HTML y datos.

`--spades-coassembly-mode global` es el valor predeterminado. `condition` requiere `--group-column`, agrupa las reads no humanas por etiqueta, ejecuta metaSPAdes por grupo y reúne contigs con procedencia. MEGAHIT conserva su coasamblea global. La validación previa rechaza grupos vacíos e identificadores de grupo que colisionen al convertirse en nombres de archivo. Este refactor no cambia esa decisión científica.

## Paralelismo y SLURM

`conf/hpc.config` conserva el executor SLURM, cuenta, partición, QoS, límites de recursos y `slurm_queue_size`. No se especifican nodos de Hypatia. Los recursos científicos posteriores se mantienen, incluido metaSPAdes con 32 CPUs y 1800 GB de RAM.

`--sra-download-max-forks N` (`sraDownloadMaxForks`, predeterminado 4) limita las tareas de adquisición concurrentes; no limita el número total de muestras. La cola y los recursos disponibles pueden reducir esa concurrencia. Los pasos FastQC, fastp y Bowtie2 reciben tareas por muestra. SRA ya no selecciona automáticamente el perfil de almacenamiento limitado.

Los perfiles `local`, `hpc`, `docker`, `conda`, `apptainer`, `singularity` y `gpu` se conservan. El modo GPU mantiene sus imágenes, argumentos y requisitos de GRES. El executor local limita la cola a una tarea cuando se activa GPU para conservar el acceso exclusivo al dispositivo.

## Reanudación y almacenamiento limitado

Añadir `--resume` al mismo comando del launcher; se traduce a `nextflow -resume`. Nextflow directo también permite `-resume <UUID>` para elegir una sesión concreta. El directorio de lanzamiento contiene la caché `.nextflow/`; `-work-dir` contiene los outputs que esa caché referencia. Ambos deben permanecer accesibles. No hay un planificador externo, una cola de pendientes ni archivos de estado por muestra.

Se conservan en `work/` los FASTQ comprimidos declarados por adquisición, reunión de runs, fastp y Bowtie2, además de los demás outputs declarados. Los reports y las reads no humanas se publican con las opciones existentes; las reads no humanas se copian por defecto. Las publicaciones son entregables y no sustituyen automáticamente a los outputs de caché. Se rechaza `publish_dir_mode=move` porque elimina archivos que Nextflow necesita para reanudar.

`--storage-constrained` selecciona `disk_efficient`: una adquisición simultánea y `storage_max_jobs=2` tareas en la cola, configurables con un archivo de configuración o parámetros nativos. Otra muestra puede estar preprocesándose mientras se adquiere el siguiente run. No se impone una barrera entre muestras o etapas.

La reducción de espacio se consigue comprimiendo cada run, evitando copias del merge de un solo run, eliminando SRA/plain FASTQ/temporales internos después de validar la salida y limitando concurrencia. **Los outputs de caché se conservan: el espacio mínimo crece con el dataset.** Este modo ya no reproduce la huella mínima del antiguo borrado de `work/`; evita que esa limpieza obligue a repetir descargas al reanudar. Revisar capacidad para reads raw, trimmed, no humanas y coasambleas antes de las 95 muestras.

`--sra-scratch-dir PATH` es opcional y usa la directiva `scratch` de Nextflow para adquisición. Nextflow crea un directorio de tarea, devuelve los outputs declarados al `work/` persistente y limpia ese scratch. En SLURM, `work/` y resultados deben ser almacenamiento compartido; comprobar montajes y capacidad del scratch en los nodos. No borrar manualmente `work/`, la caché, ni outputs publicados que se quieran conservar mientras se necesite reanudar.

Las sesiones de la arquitectura anterior tienen nombres de procesos y límites de caché diferentes; no se garantiza reutilización de su caché. Los resultados o FASTQ antiguos no se importan basándose sólo en su nombre. Usar una ubicación de resultados nueva para la primera validación del refactor. El refactor local no modifica ningún archivo de Hypatia.

## Configuración y recursos externos

La configuración de bases de datos existente debe definir `host_bowtie2_index`, `checkm2_db`, `gunc_db`, `gtdbtk_db`, `phylophlan_db`, `genemark_home`, `genemark_key`, `eggnog_db` e `interproscan_data`. Se conserva `assets/phylophlan_iqtree.cfg` y los ajustes específicos ya existentes de InterProScan. La instalación/licencia de GeneMark es externa.

`bin/prepare_databases.sh` y `bin/prepare_containers.sh` son utilidades operativas conservadas para preparar recursos por separado. Consultar su `--help`. No forman una cola de muestras y no se ejecutan automáticamente para resolver inputs. La adquisición utiliza SRA Toolkit, pigz, Bash y AWK. La receta SRA añade gawk explícitamente: verificar que una imagen previamente cacheada incluya `awk` y reconstruirla si falta. No se construyeron ni descargaron imágenes durante el refactor.

El launcher no necesita Python en el host para ejecutar Nextflow. Conserva el bloqueo del directorio de resultados, los códigos de salida y los informes nativos por invocación en `pipeline_info/invocations/`. `reportsDir` controla el destino de trace/report/timeline/DAG nativos; la configuración HPC mantiene trace/report/timeline desactivados. Las métricas GPU siguen publicándose en `pipeline_info/resources/gpu_tasks/` y conservan `--gpu-telemetry-interval`. Las utilidades de preparación de recursos mantienen sus propios requisitos de Python.

## Resultados y estructura

```text
assets/raw/         manifests de 8 y 95 BioSamples
bin/                transformaciones científicas y utilidades operativas
conf/               perfiles, recursos y publicación
containers/         recetas de imágenes mantenidas por el proyecto
lib/                validación de metadata y agrupación Groovy
modules/core/       herramientas bioinformáticas reutilizables
modules/local/      integración, catálogos, metadata SRA y reunión de runs
subworkflows/local/ entrada, preprocesamiento y etapas científicas
workflows/          selección de entrada local/SRA
```

Los directorios de resultados `01_quality_control_and_filtering` a `06_global_processing_evaluation` conservan su estructura. `pipeline_info/` añade la tabla común `sample_metadata.tsv`, `sra_exclusions.tsv`, el RunInfo recuperado y logs SRA; conserva versiones de software e informes de ejecución. Las columnas adicionales del manifest siguen presentes en la tabla de metadata. La tabla larga de abundancia mantiene el contrato científico previo, incluidos `sample_id` y `group`.

Los siete auxiliares Python restantes y la sustitución del filtro FASTA por AWK se documentan en [docs/python_reduction_review.md](docs/python_reduction_review.md). El primer refactor se describe en [docs/refactor_review.md](docs/refactor_review.md). La arquitectura detallada está en [docs/architecture.md](docs/architecture.md); las auditorías de recursos, compresión y capacidades GPU conservan sus referencias científicas.

## Validación ligera

Sin ejecutar herramientas ni descargar reads o bases de datos:

```bash
nextflow -version
nextflow config -profile hpc,apptainer -flat
nextflow config -profile hpc,apptainer,disk_efficient -flat
nextflow inspect . -profile hpc,apptainer
bash -n metagenomics_pipeline.sh
git diff --check
```

`nextflow inspect` no debe llevar `-concretize`, que construiría imágenes. No se ofrece un `-stub-run` integral: los módulos científicos posteriores no tienen stubs completos. La primera ejecución científica será la validación de ocho muestras en Hypatia, después de revisar el diff.
