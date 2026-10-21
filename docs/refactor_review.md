# Revisión del refactor SRA nativo de Nextflow

Primer refactor: 2026-10-20, commit `2a6b2e4`. Actualizado el 2026-10-21 tras la segunda revisión de Python, sobre `refactor/native-nextflow-sra`. Los cambios de esta segunda revisión quedan sin commit, sin staging y sin push. No se ejecutó el pipeline científico, no se descargaron reads, bases de datos o imágenes, y no se modificó Hypatia.

## Auditoría y decisión de arquitectura

Se inspeccionaron el entrypoint, workflow dispatcher, todos los subworkflows y contratos de módulos, configuración de recursos/publicación/perfiles, scripts de `bin/`, entornos/recetas, manifests, esquema, documentación, CI y colección de pruebas. No existían carpetas productivas `scripts/` ni `launchers/`: el launcher era el archivo de la raíz. Se buscaron referencias antes de retirar auxiliares.

| Familia revisada | Papel comprobado | Decisión |
| --- | --- | --- |
| `main.nf`, `workflows/` | Entrada única y selección de modalidad | Conservar entrypoint; sustituir dispatch de etapas internas por entrada local/SRA |
| `subworkflows/local/sra*` anteriores | Descubrimiento, reconciliación, una muestra por sesión, análisis global separado | Sustituir por entrada y adquisición DSL2 con channels del conjunto completo |
| Módulos locales de input/estado anteriores | Conversión de tablas, commit de FASTQ y JSON de finalización | Retirar cuando sus consumidores fueron sustituidos |
| `modules/core/sratools` | Invocación de un controlador Python por muestra | Tarea por run con herramientas directas y validación FASTQ AWK |
| Preprocesamiento | FastQC/fastp/FastQC/Bowtie2 y joins por ID | Sin cambios científicos ni cambios en el propio subworkflow |
| Construcción MAGs y análisis posteriores | Coasambleas, binning, refinamiento, taxonomía, anotación, abundancia, MultiQC | Conservar módulos/subworkflows científicos y parámetros |
| `conf/` | Recursos y publicación, más serialización global de disco | Retirar selectores de procesos eliminados; limitar adquisición/cola sólo si se solicita |
| `tests/` | Unit tests, contratos estáticos, renderizados/mocks de regresiones del desarrollo y fixtures | Ningún consumidor productivo; retirar la colección según la solicitud |
| `bin/prepare_databases.sh` | Preparación/validación de recursos externos y configuración | Conservar; no ejecutar |
| `bin/prepare_containers.sh` | Inventario/preparación/validación de imágenes | Conservar; no ejecutar |
| `.github/workflows/build-containers.yml` | Publicación de imágenes | Conservar sin cambios; no disparar |
| Auditorías de recursos/compresión/GPU | Evidencia científica y operativa | Conservar; actualizar únicamente referencias de la arquitectura SRA |
| Artefactos locales ignorados | Bytecode, logs y reports sin tareas | Limpiar artefactos identificados; conservar recursos/entornos y cachés preexistentes ajenos a la validación |

La arquitectura anterior no se retiró a ciegas: el launcher ejecutaba descubrimiento, reconciliación inicial, un `nextflow run` por muestra pendiente, reconciliación final y una sesión global. `acquire_sra_sample.py` secuenciaba las herramientas; `manage_sra_checkpoints.py` validaba, copiaba y reconciliaba pares no humanos, mantenía JSON/complete y sellaba outputs. El launcher borraba el trabajo de cada muestra y luego el trabajo global. Por eso la caché nativa sola no era suficiente: sus archivos habían desaparecido. Casi toda ejecución SRA activaba además una cola global de una tarea.

La arquitectura nueva lee y valida el TSV completo con Nextflow/Groovy, resuelve metadata una vez y emite todas las muestras. `SRA_ACQUISITION` expande sus runs a tareas independientes, reúne cada BioSample con `groupKey`/`groupTuple` y emite `tuple(meta, reads)` al preprocesamiento común. Un run tardío de otra muestra no bloquea una muestra lista. La reunión de reads de la coasamblea posterior sigue siendo una dependencia científica necesaria.

No se concatenan muestras para Bowtie2. `meta` conserva BioSample, grupo, proyecto, todos sus runs y las columnas originales del input. `PD` y `CT` permanecen sin traducciones; otros nombres de grupos son válidos. El subset tiene 4 PD y 4 CT; el completo, confirmado como 95 por el usuario, tiene 46 PD y 49 CT. El subset fue renombrado conservando exactamente su contenido; el manifest completo no se modificó.

## Python conservado y consumidores actuales

| Archivo | Consumidor real | Justificación funcional |
| --- | --- | --- |
| `bin/fetch_sra_runinfo.py` | `FETCH_SRA_RUNINFO` | Adaptador HTTP NCBI ESearch history/EFetch, reintentos de red y materialización atómica de RunInfo bruto. No lee el TSV ni selecciona, agrupa, programa o controla muestras. Conserva el acceso a metadata sin introducir EDirect/curl+jq como nuevas dependencias. |
| `bin/select_high_quality_mags.py` | `SELECT_HIGH_QUALITY_MAGS` | Cruza bins con CheckM2 y aplica los umbrales científicos estrictos; escribe tablas para dRep. |
| `bin/combine_mag_catalogs.py` | `COMBINE_MAG_CATALOGS` | Identificadores únicos, calidad y procedencia al integrar catálogos de dos ensambladores. |
| `bin/finalize_mag_catalog.py` | `FINALIZE_MAG_CATALOG` | Materializa representantes y comprueba correspondencia de FASTA/calidad/procedencia. |
| `bin/spades_coassembly.py` | `MERGE_SPADES_ASSEMBLIES` | Sólo merge científico de contigs por condición y procedencia. Se eliminaron validación de inputs y estado de estrategia externos. |
| `bin/integrate_functional_annotations.py` | `INTEGRATE_ANNOTATIONS` | Integra FASTA/GFF/eggNOG/InterPro por gen y MAG; transformación científica multiformato. |
| `bin/normalize_coverm_abundance.py` | `NORMALIZE_ABUNDANCE` | Valida métricas y transforma abundancia a tabla larga con grupos, manteniendo valores científicos. |


El launcher no requiere Python para ejecutar el workflow. Los informes nativos y las métricas GPU se conservan. La preparación de recursos mantiene sus dependencias independientes. Ningún auxiliar mantiene una cola de muestras.

## Parámetros retirados

Del contrato Nextflow: `sraManifest`, `sraCheckpointManifest`, `sraSampleMetadata`, `sraSampleId`, `sraCheckpointDir`, `sraCacheDir`, `sraTempDir`, `sraStateDir`, `sraRequireComplete`, `executionStage`.

Del launcher: `--sra-checkpoint-dir`, `--keep-sra-checkpoints`, `--sra-cache-dir`, `--sra-temp-dir`. Ya no se generan registros externos de sesiones de reanudación ni marcadores de éxito por muestra/proyecto que actúen como controladores. Las menciones a nombres retirados en este informe son inventario histórico, no referencias ejecutables.

## Parámetros nuevos o modificados

| Parámetro | Comportamiento |
| --- | --- |
| `sraRunInfo` / `--sra-run-info` | CSV NCBI existente opcional; evita consulta online de metadata, manteniendo validación científica |
| `sraDownloadMaxForks` / `--sra-download-max-forks` | Máximo de adquisiciones concurrentes, predeterminado 4 |
| `storage_max_jobs` | Cola del perfil `disk_efficient`, predeterminado 2 |
| `sraScratchDir` / `--sra-scratch-dir` | Opcional; directiva nativa `scratch` de Nextflow, ya no raíz de un controlador de muestras |
| `sraContainerOptions` | Opciones de runtime opcionales; no genera binds para estado/cache/temp externos |
| `--storage-constrained` | Perfil explícito: una adquisición y cola configurable; SRA deja de activarlo automáticamente |
| `--resume` | Una invocación con `-resume` nativo, sin reconciliación o registro paralelo de UUIDs |
| `--work-dir` | Trabajo persistente común, independiente del scratch; no se borra al terminar muestras |
| `groupColumn` | Mantiene columna configurable; `group` se detecta si existe y no se indica otra |
| `publish_dir_mode` | Se rechaza `move` para evitar eliminar outputs necesarios para caché/resume |

Se conservan `sraProject`, `sraSamples`, plataformas, máximo de prefetch, contacto NCBI, imágenes científicas, parámetros de grupos/coasamblea, recursos SLURM y opciones científicas de las etapas posteriores.

## Reanudación y almacenamiento

`-resume` reutiliza tareas cuyos inputs, metadata, código y outputs cacheados siguen siendo válidos. Se conservan `.nextflow/` y `work/`; la posición en el manifest no se introduce en el hash por muestra. Los FASTQ raw/merged/trimmed/no humanos declarados no se borran desde otros procesos. Los FASTQ publicados son entregables explícitos y no se adoptan como caché sólo por su nombre.

En cada adquisición se eliminan únicamente SRA, FASTQ sin comprimir y temporales de esa tarea después de validar mates y gzip. Para un BioSample de un run, el merge usa enlaces; varios runs se reúnen por mate en orden determinista. El scratch nativo devuelve outputs al trabajo persistente antes de limpiarse.

El modo limitado reduce temporales simultáneos, con una adquisición y dos tareas en cola por defecto. Permite intercalar adquisición/preprocesamiento; no espera a completar toda una muestra para admitir las demás. **Conserva más datos que el borrado anterior:** requiere espacio para outputs cacheados de todas las muestras y coasambleas. No ofrece una garantía de espacio independiente del número de muestras. Esta elección evita rehacer descargas al reanudar tareas válidas.

## Comandos para Hypatia (no ejecutados)

Definir rutas reales de la configuración ya preparada y del almacenamiento persistente:

```bash
DB_CONFIG=/ruta/real/metagenomics_databases.config
RUN_ROOT=/ruta/persistente/PRJNA782492
```

Ocho muestras:

```bash
bash metagenomics_pipeline.sh --hpc --apptainer --run \
  --sra-project PRJNA782492 \
  --sra-samples assets/raw/PRJNA782492_subset8.tsv \
  --group-column group --database-config "$DB_CONFIG" \
  --outdir "$RUN_ROOT/results" --work-dir "$RUN_ROOT/work"
```

Conjunto completo, 95 muestras:

```bash
bash metagenomics_pipeline.sh --hpc --apptainer --run \
  --sra-project PRJNA782492 \
  --sra-samples assets/raw/PRJNA782492_metagenomic_samples.tsv \
  --group-column group --database-config "$DB_CONFIG" \
  --outdir "$RUN_ROOT/results" --work-dir "$RUN_ROOT/work"
```

Añadir `--resume` al reanudar o ampliar una ejecución propia del refactor. Añadir `--storage-constrained` si se desea el perfil limitado. Para comparar corridas independientes utilizar otro `RUN_ROOT`; el código no contiene condiciones por cantidad de muestras.

## Validación del primer refactor y límites

Las comprobaciones siguientes corresponden al primer refactor, commit `2a6b2e4`, con Nextflow 26.04.6 en WSL. La validación de los cambios actuales se registra en [python_reduction_review.md](python_reduction_review.md). No se usó `-stub-run` integral ni se ejecutaron módulos científicos reales.

| Comprobación | Resultado |
| --- | --- |
| `nextflow -version` | 26.04.6 |
| `nextflow config -flat` | Perfiles local/docker, hpc/apptainer, hpc/apptainer/disk_efficient y hpc/apptainer/gpu válidos |
| Configuración de concurrencia | HPC: cola 100 y adquisiciones 4; disk: cola 2 y adquisiciones 1 |
| `nextflow inspect . -profile hpc,apptainer` sin concretize | 33 definiciones; procesos nuevos presentes y retirados ausentes |
| Metadata/channels aislados | TSV 8/95, grupos 4+4 del subset, grupos genéricos, columnas extra, duplicados/campos inválidos, selección WGS y orden de runs correctos |
| Metadata al ampliar manifest | Identidad y metadata por muestra estables entre subset y conjunto completo |
| Subworkflow SRA real con ejecutables simulados | Dos BioSamples, tres runs; FASTQ de cuatro bases, todos los mates correctos, grupo conservado, scratch limpio |
| Repetición simulada con `-resume` | Cinco tareas CACHED, cero adquisiciones adicionales |
| Ampliación simulada a tres BioSamples/cuatro runs | Cinco tareas CACHED y sólo dos tareas nuevas; una adquisición adicional |
| AWK | Ocho escenarios válidos/erróneos: estructura, IDs, orientación de mates, longitud, truncamiento y número de registros |
| Adaptador HTTP NCBI | Respuestas en memoria simuladas; cero conexiones de red; CSV quoted y materialización correctos |
| Launcher | Dry-run único para 8/95/local, perfiles, resume y argumentos con `=`; sin escrituras |
| Sintaxis | Python AST, Bash `-n`, JSON y equivalencia del esquema con parámetros declarados |
| Referencias | Includes y scripts existentes, sin consumidores productivos de auxiliares/etapas retiradas |
| Protección científica | Módulos core distintos de SRA y subworkflows posteriores/preprocesamiento iguales a HEAD; MultiQC conserva `multiqc_output/multiqc_report_data` |
| Git | `git diff --check` sin errores; rama actual conservada; índice sin cambios |

Las pruebas sintéticas y sus herramientas simuladas son artefactos temporales locales, no una nueva colección de tests productivos. Su evidencia se resume aquí y se conserva temporalmente fuera del repositorio durante la limpieza.

Antes de ejecutar en Hypatia revisar: capacidad para outputs cacheados; acceso compartido a `work/` y resultados; montajes/capacidad del scratch si se usa; rutas reales de bases de datos/licencia e índices; cuenta/partición/QoS de SLURM; disponibilidad de las imágenes y del ejecutable `awk` en la imagen SRA. Las recetas Conda/imagen ahora incluyen gawk explícitamente; una imagen antigua que no tenga AWK debe reconstruirse. No se descargaron/construyeron imágenes para verificar su contenido aquí.

No se garantiza reanudar caché de la arquitectura anterior: cambiaron procesos y límites de sesiones, y el diseño anterior eliminaba `work/`. No se importaron FASTQ antiguos por nombre ni se modificaron resultados remotos. La primera validación real debe usar las ocho muestras con una ubicación nueva. Se mantiene la desactivación previa de trace/report/timeline en HPC, con posibles métricas no disponibles. La equivalencia científica se sustenta en preservar código/argumentos y contratos; la ejecución biológica real está pendiente por instrucción del usuario.


El inventario y los snapshots Git del primer refactor pueden consultarse en el commit `2a6b2e4`. Esta documentación describe el árbol actual; la segunda revisión se documenta en [python_reduction_review.md](python_reduction_review.md).
