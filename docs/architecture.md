# Arquitectura del pipeline

## Contratos de entrada

`main.nf` invoca `METAGENOMICS`, que selecciona `LOCAL_INPUT` o `SRA_INPUT`. Ambos generan `tuple(meta, [read_1, read_2])`, ejecutan `QUALITY_CONTROL_AND_FILTERING` y conectan sus outputs a `METAGENOMICS_GLOBAL`.

`InputMetadata` valida el conjunto completo antes de emitir muestras: cabeceras, número de columnas, identidades, grupos y colisiones de identificadores de coasamblea. Conserva todas las columnas originales bajo `meta.metadata`. `group` es la categoría normalizada por trim, sin una taxonomía fija de etiquetas. `sample_metadata.tsv` contiene las columnas comunes `sample_id`, `biosample_accession`, `group`, `run_accessions`, `project_accession`, `sample_order` y las columnas adicionales del input.

La posición de una fila sólo pertenece a la tabla del conjunto; no altera el hash de tareas por muestra cuando se amplía una selección. Las rutas FASTQ relativas se resuelven respecto al CSV original.

## SRA nativo de Nextflow

1. `splitCsv(header:false, sep:'\t')` lee todo el TSV. Una reunión inicial de metadata permite rechazar un conjunto inválido antes de descargar reads. No es una barrera entre etapas científicas.
2. `FETCH_SRA_RUNINFO` obtiene metadata bruta mediante NCBI ESearch history/EFetch. Es una tarea cacheable por BioProject. `sraRunInfo` permite usar un CSV existente y omitir esa consulta.
3. `SraMetadata` selecciona los BioSamples explícitos, valida elegibilidad y conserva todos sus runs aptos. Rechaza runs duplicados y BioSamples sin datos elegibles. Conserva los filtros previos de acceso público, BioProject, PAIRED/WGS/METAGENOMIC, plataformas y contadores.
4. El channel completo de muestras se expande a tareas `SRATOOLS_ACQUIRE(meta, run)`. Cada tarea usa `prefetch`, `vdb-validate`, `fasterq-dump --split-files --skip-technical`, validación AWK de ambos mates, pigz y CRC. Las estructuras SRA y FASTQ sin comprimir sólo viven en el directorio propio de la tarea.
5. `groupKey(meta, run_count)` y `groupTuple` liberan cada muestra al completar sus propios runs. Los pares se ordenan por accesión y por mate, independientemente del orden de finalización.
6. `MERGE_SRA_RUNS` combina únicamente runs de ese BioSample. Un solo run usa enlaces simbólicos; varios runs se descomprimen y recomprimen por mate a un flujo gzip. No se combinan muestras distintas en esta etapa.
7. Los pares con su metadata alimentan el preprocesamiento común sin otra invocación de Nextflow.

El único Python de entrada SRA es el adaptador HTTP de metadata. No lee el TSV de selección, no filtra muestras, no ejecuta herramientas, no maneja estados y no programa tareas. El scheduler de Nextflow posee todas las dependencias y reintentos de procesos.

## Preprocesamiento compartido

FastQC raw -> fastp -> FastQC clean -> Bowtie2 host removal conserva herramientas, versiones y argumentos científicos. Los joins por `meta.id` hacen que cada muestra espere únicamente su propio FastQC. Bowtie2 recibe un par por muestra y conserva `meta` en las reads no humanas y logs. No existe una reunión global de FASTQ antes de Bowtie2.

## Análisis posterior conservado

`METAGENOMICS_GLOBAL` recibe reads no humanas, reports, versiones y metadata. Conserva los subworkflows de construcción de MAGs, taxonomía/filogenómica, anotación funcional, abundancia y evaluación global.

- MEGAHIT conserva coasamblea global. metaSPAdes conserva modos global y por condición. Sus reuniones de muestras son necesarias para la coasamblea. El merge científico de contigs por condición preserva grupos e identificadores originales.
- Cada rama conserva MetaQUAST, filtrado de contigs, CoverM, cuatro binners independientes (COMEBin, MetaBAT2, SemiBin2, Vamb), DAS Tool, CheckM2/GUNC, selección HQ y dRep 99%/95%.
- Se combinan los catálogos con identificadores y procedencia; dRep selecciona el catálogo final, que vuelve a evaluarse.
- GTDB-Tk y PhyloPhlAn/IQ-TREE usan ese catálogo. GeneMarkS-2 produce proteínas/CDS/GFF; eggNOG-mapper e InterProScan se integran por MAG.
- CoverM usa la colección de muestras sin perder sus identificadores; la transformación a tabla larga mantiene la asignación a grupos del contrato previo.
- MultiQC reúne reports y conserva la ruta `multiqc_output/multiqc_report_data` requerida por 1.35.

Los módulos científicos posteriores, sus argumentos, versiones y reglas de publicación no se reescribieron. Se mantiene la limpieza interna ya existente de árboles nativos temporales después de materializar los outputs declarados.

## Caché y almacenamiento

La reanudación usa únicamente caché nativa y `-resume`. El launcher lanza una sola ejecución y no elimina outputs de otros procesos. Las tareas SRA eliminan sus propios temporales después de validar los outputs; Nextflow conserva los FASTQ declarados y sus dependencias en `work/`.

La directiva `scratch` opcional de adquisición copia los outputs declarados a `work/` persistente antes de limpiar el scratch. No sustituye a `work/`. Los entregables publicados son copias/enlaces configurados; no hay descubrimiento automático de FASTQ publicados para decidir qué tareas omitir.

El perfil `disk_efficient` limita a una adquisición y a `storage_max_jobs` tareas en cola (predeterminado 2). El resto del dataflow conserva sus dependencias por muestra y el executor puede intercalar etapas. SRA no activa el perfil por defecto.

Esta política preserva reanudación sin descargas repetidas para tareas cacheadas válidas, a cambio de conservar raw/trimmed/non-host declarados. El consumo mínimo aumenta con el número de muestras. Limpiar `work/` manualmente invalida esa garantía; no se presenta el perfil como una garantía de capacidad total.

## Configuración, SLURM y observación

`nextflow.config` incluye base, recursos y publicación, y después los perfiles elegidos. Se conservan recursos y parámetros SLURM, sin nodos fijados. `sraDownloadMaxForks` limita adquisiciones; `slurm_queue_size` y `local_max_jobs` controlan el executor normal. `disk_efficient` es explícito.

El modo GPU mantiene imágenes, versiones y GRES existentes. La cola local de una tarea al habilitar GPU conserva el acceso exclusivo al dispositivo. Los límites de concurrencia no crean dependencias científicas artificiales.

Los tres auxiliares opcionales de recursos observan almacenamiento, `sacct` y trazas para informes. Sus datos no se usan para seleccionar, omitir o reordenar muestras. No hay registros externos de reanudación. La configuración HPC previa desactiva trace/report/timeline por compatibilidad con contenedores mínimos, y puede dejar métricas sin datos; esta limitación sigue siendo explícita.

Las utilidades de bases de datos y contenedores se conservan por ser operativas. No se ejecutan durante la lectura del TSV ni durante las verificaciones estáticas.

## Migración y validación

La primera validación real utiliza `PRJNA782492_subset8.tsv`. La ejecución completa utiliza `PRJNA782492_metagenomic_samples.tsv`, con 95 BioSamples confirmados por el usuario. No hay condiciones por cantidad en el código productivo.

El diseño anterior repartía una ejecución científica entre varias sesiones y borraba directorios de trabajo después de cada muestra. El diseño nuevo cambia esos límites de caché. No se garantiza compatibilidad de reanudación con sesiones antiguas, ni se importan resultados sólo por coincidencia de nombres. Se recomienda una ubicación nueva para la primera validación.

La validación local se limita a inspección/configuración de Nextflow, checks de metadata y channels sintéticos, shell, referencias y diff. El informe [refactor_review.md](refactor_review.md) registra resultados y archivos afectados. No se ejecuta el pipeline científico ni se modifican recursos remotos.
