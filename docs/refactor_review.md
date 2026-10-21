# Revisión del refactor SRA nativo de Nextflow

Fecha: 2026-10-20. Rama local: `refactor/native-nextflow-sra`. Cambios sin commit, sin staging y sin push. No se ejecutó el pipeline científico, no se descargaron reads, bases de datos o imágenes, y no se modificó Hypatia.

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

## Python conservado y consumidores

| Archivo | Consumidor real | Justificación funcional |
| --- | --- | --- |
| `bin/fetch_sra_runinfo.py` | `FETCH_SRA_RUNINFO` | Adaptador HTTP NCBI ESearch history/EFetch, reintentos de red y materialización atómica de RunInfo bruto. No lee el TSV ni selecciona, agrupa, programa o controla muestras. Conserva el acceso a metadata sin introducir EDirect/curl+jq como nuevas dependencias. |
| `bin/filter_fasta_by_length.py` | `FILTER_CONTIGS` en binning | Transformación científica streaming de FASTA y estadísticas; conserva umbral y formato previos. |
| `bin/select_high_quality_mags.py` | `SELECT_HIGH_QUALITY_MAGS` | Cruza bins con CheckM2 y aplica los umbrales científicos estrictos; escribe tablas para dRep. |
| `bin/combine_mag_catalogs.py` | `COMBINE_MAG_CATALOGS` | Identificadores únicos, calidad y procedencia al integrar catálogos de dos ensambladores. |
| `bin/finalize_mag_catalog.py` | `FINALIZE_MAG_CATALOG` | Materializa representantes y comprueba correspondencia de FASTA/calidad/procedencia. |
| `bin/spades_coassembly.py` | `MERGE_SPADES_ASSEMBLIES` | Sólo merge científico de contigs por condición y procedencia. Se eliminaron validación de inputs y estado de estrategia externos. |
| `bin/integrate_functional_annotations.py` | `INTEGRATE_ANNOTATIONS` | Integra FASTA/GFF/eggNOG/InterPro por gen y MAG; transformación científica multiformato. |
| `bin/normalize_coverm_abundance.py` | `NORMALIZE_ABUNDANCE` | Valida métricas y transforma abundancia a tabla larga con grupos, manteniendo valores científicos. |
| `bin/monitor_storage.py` | Observador opcional del launcher | Mide almacenamiento sin seguir enlaces y registra picos; nunca decide qué se procesa. Se retiraron categorías y handshake de limpieza de la arquitectura anterior. |
| `bin/collect_slurm_accounting.py` | Observador opcional del launcher | Consulta `sacct` para informes; no envía ni administra jobs científicos. |
| `bin/summarize_resources.py` | Informe opcional del launcher | Agrega trazas/mediciones. Se adaptaron scopes y listas de runs; no controla el workflow. |

Los ocho auxiliares científicos/de metadata evitan reescribir etapas posteriores fuera del alcance. Los tres observadores se desactivan con `--no-resource-telemetry`, sin Python obligatorio en el host de lanzamiento. El JSON de resultado usado por esos informes es observacional: nunca autoriza omitir tareas o limpiar outputs. El Python embebido conservado en los módulos GPU comprueba CUDA/PyTorch; el de preparación de recursos maneja sus inventarios/configuración. Ninguno mantiene una cola de muestras.

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
| `--no-resource-telemetry` | Desactiva los tres observadores opcionales del launcher |
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

## Validación realizada y límites

Todas las comprobaciones finales pasaron con Nextflow 26.04.6 en WSL. No se usó `-stub-run` integral ni se ejecutaron módulos científicos reales.

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
| Observadores | Muestreo único y resumen con trazas sintéticas y tabla de metadata nueva |
| Sintaxis | Python AST, Bash `-n`, JSON y equivalencia del esquema con parámetros declarados |
| Referencias | Includes y scripts existentes, sin consumidores productivos de auxiliares/etapas retiradas |
| Protección científica | Módulos core distintos de SRA y subworkflows posteriores/preprocesamiento iguales a HEAD; MultiQC conserva `multiqc_output/multiqc_report_data` |
| Git | `git diff --check` sin errores; rama actual conservada; índice sin cambios |

Las pruebas sintéticas y sus herramientas simuladas son artefactos temporales locales, no una nueva colección de tests productivos. Su evidencia se resume aquí y se conserva temporalmente fuera del repositorio durante la limpieza.

Antes de ejecutar en Hypatia revisar: capacidad para outputs cacheados; acceso compartido a `work/` y resultados; montajes/capacidad del scratch si se usa; rutas reales de bases de datos/licencia e índices; cuenta/partición/QoS de SLURM; disponibilidad de las imágenes y del ejecutable `awk` en la imagen SRA. Las recetas Conda/imagen ahora incluyen gawk explícitamente; una imagen antigua que no tenga AWK debe reconstruirse. No se descargaron/construyeron imágenes para verificar su contenido aquí.

No se garantiza reanudar caché de la arquitectura anterior: cambiaron procesos y límites de sesiones, y el diseño anterior eliminaba `work/`. No se importaron FASTQ antiguos por nombre ni se modificaron resultados remotos. La primera validación real debe usar las ocho muestras con una ubicación nueva. Se mantiene la desactivación previa de trace/report/timeline en HPC, con posibles métricas no disponibles. La equivalencia científica se sustenta en preservar código/argumentos y contratos; la ejecución biológica real está pendiente por instrucción del usuario.

## Inventario final de archivos

Retirados: 38; renombrado: 1; nuevos incluyendo el manifest renombrado: 10; modificados: 23. Python retirados: 16; Python productivos restantes: 11.

### Archivos retirados

- `bin/acquire_sra_sample.py`
- `bin/check_samplesheet.py`
- `bin/manage_sra_checkpoints.py`
- `bin/resolve_sra_project.py`
- `modules/local/check_samplesheet/environment.yml`
- `modules/local/check_samplesheet/main.nf`
- `modules/local/sra_preprocessing/environment.yml`
- `modules/local/sra_preprocessing/main.nf`
- `modules/local/sra_project/environment.yml`
- `modules/local/sra_project/main.nf`
- `subworkflows/local/sra_checkpoint_reconciliation/main.nf`
- `subworkflows/local/sra_global/main.nf`
- `subworkflows/local/sra_project_discovery/main.nf`
- `subworkflows/local/sra_sample_preprocessing/main.nf`
- `tests/fixtures/resource_accounting/human_trace.tsv`
- `tests/fixtures/resource_accounting/raw_trace.tsv`
- `tests/fixtures/resource_accounting/sra_manifest.tsv`
- `tests/fixtures/resource_accounting/storage_usage_timeseries.tsv`
- `tests/fixtures/resource_accounting/task_workdir_peaks.tsv`
- `tests/fixtures/resource_accounting/task_workdir_timeseries.tsv`
- `tests/fixtures/resource_accounting/trace_registry.tsv`
- `tests/fixtures/runinfo_invalid.csv`
- `tests/fixtures/runinfo_selection.csv`
- `tests/fixtures/runinfo_valid.csv`
- `tests/fixtures/sra_samples_invalid.tsv`
- `tests/fixtures/sra_samples_valid.tsv`
- `tests/test_check_samplesheet.py`
- `tests/test_coverm_contig_awk.py`
- `tests/test_coverm_mag_abundance_awk.py`
- `tests/test_genemarks2_rendering.py`
- `tests/test_launcher_lock_contract.py`
- `tests/test_monitor_storage.py`
- `tests/test_normalize_abundance.py`
- `tests/test_phylophlan_rendering.py`
- `tests/test_slurm_accounting.py`
- `tests/test_sra_helpers.py`
- `tests/test_static_contracts.py`
- `tests/test_summarize_resources.py`

### Renombrado

`assets/raw/test_metagenomic_samples.tsv` -> `assets/raw/PRJNA782492_subset8.tsv`. Mismo hash Git; ocho muestras sin cambios.

### Archivos nuevos (incluido el renombrado)

- `assets/raw/PRJNA782492_subset8.tsv`
- `bin/fetch_sra_runinfo.py`
- `bin/validate_fastq_pair.awk`
- `docs/refactor_review.md`
- `lib/InputMetadata.groovy`
- `lib/SraMetadata.groovy`
- `modules/local/merge_sra_runs/main.nf`
- `modules/local/sra_runinfo/environment.yml`
- `modules/local/sra_runinfo/main.nf`
- `subworkflows/local/sra_acquisition/main.nf`

### Archivos modificados

- `.gitignore`
- `README.md`
- `bin/monitor_storage.py`
- `bin/spades_coassembly.py`
- `bin/summarize_resources.py`
- `conf/disk_efficient.config`
- `conf/modules.config`
- `conf/resources.config`
- `containers/sra-tools/Dockerfile`
- `containers/sra-tools/README.md`
- `containers/sra-tools/environment.yml`
- `docs/architecture.md`
- `docs/compression_audit.tsv`
- `docs/original_resource_mapping.tsv`
- `metagenomics_pipeline.sh`
- `modules/core/sratools/environment.yml`
- `modules/core/sratools/main.nf`
- `modules/core/sratools/meta.yml`
- `nextflow.config`
- `nextflow_schema.json`
- `subworkflows/local/local_input/main.nf`
- `subworkflows/local/sra_input/main.nf`
- `workflows/metagenomics.nf`

### Todos los Python retirados

- `bin/acquire_sra_sample.py`
- `bin/check_samplesheet.py`
- `bin/manage_sra_checkpoints.py`
- `bin/resolve_sra_project.py`
- `tests/test_check_samplesheet.py`
- `tests/test_coverm_contig_awk.py`
- `tests/test_coverm_mag_abundance_awk.py`
- `tests/test_genemarks2_rendering.py`
- `tests/test_launcher_lock_contract.py`
- `tests/test_monitor_storage.py`
- `tests/test_normalize_abundance.py`
- `tests/test_phylophlan_rendering.py`
- `tests/test_slurm_accounting.py`
- `tests/test_sra_helpers.py`
- `tests/test_static_contracts.py`
- `tests/test_summarize_resources.py`

Los cuatro auxiliares productivos retirados implementaban adquisición/control externo, reconciliación o adaptación de inputs ahora cubierta por Nextflow/Groovy. Los otros doce eran exclusivamente la colección de pruebas de desarrollo. Ningún archivo Python retirado era necesario para preparar bases de datos o para las transformaciones científicas posteriores conservadas.

## Estado Git al terminar

```text
 M .gitignore
 M README.md
 D assets/raw/test_metagenomic_samples.tsv
 D bin/acquire_sra_sample.py
 D bin/check_samplesheet.py
 D bin/manage_sra_checkpoints.py
 M bin/monitor_storage.py
 D bin/resolve_sra_project.py
 M bin/spades_coassembly.py
 M bin/summarize_resources.py
 M conf/disk_efficient.config
 M conf/modules.config
 M conf/resources.config
 M containers/sra-tools/Dockerfile
 M containers/sra-tools/README.md
 M containers/sra-tools/environment.yml
 M docs/architecture.md
 M docs/compression_audit.tsv
 M docs/original_resource_mapping.tsv
 M metagenomics_pipeline.sh
 M modules/core/sratools/environment.yml
 M modules/core/sratools/main.nf
 M modules/core/sratools/meta.yml
 D modules/local/check_samplesheet/environment.yml
 D modules/local/check_samplesheet/main.nf
 D modules/local/sra_preprocessing/environment.yml
 D modules/local/sra_preprocessing/main.nf
 D modules/local/sra_project/environment.yml
 D modules/local/sra_project/main.nf
 M nextflow.config
 M nextflow_schema.json
 M subworkflows/local/local_input/main.nf
 D subworkflows/local/sra_checkpoint_reconciliation/main.nf
 D subworkflows/local/sra_global/main.nf
 M subworkflows/local/sra_input/main.nf
 D subworkflows/local/sra_project_discovery/main.nf
 D subworkflows/local/sra_sample_preprocessing/main.nf
 D tests/fixtures/resource_accounting/human_trace.tsv
 D tests/fixtures/resource_accounting/raw_trace.tsv
 D tests/fixtures/resource_accounting/sra_manifest.tsv
 D tests/fixtures/resource_accounting/storage_usage_timeseries.tsv
 D tests/fixtures/resource_accounting/task_workdir_peaks.tsv
 D tests/fixtures/resource_accounting/task_workdir_timeseries.tsv
 D tests/fixtures/resource_accounting/trace_registry.tsv
 D tests/fixtures/runinfo_invalid.csv
 D tests/fixtures/runinfo_selection.csv
 D tests/fixtures/runinfo_valid.csv
 D tests/fixtures/sra_samples_invalid.tsv
 D tests/fixtures/sra_samples_valid.tsv
 D tests/test_check_samplesheet.py
 D tests/test_coverm_contig_awk.py
 D tests/test_coverm_mag_abundance_awk.py
 D tests/test_genemarks2_rendering.py
 D tests/test_launcher_lock_contract.py
 D tests/test_monitor_storage.py
 D tests/test_normalize_abundance.py
 D tests/test_phylophlan_rendering.py
 D tests/test_slurm_accounting.py
 D tests/test_sra_helpers.py
 D tests/test_static_contracts.py
 D tests/test_summarize_resources.py
 M workflows/metagenomics.nf
?? assets/raw/PRJNA782492_subset8.tsv
?? bin/fetch_sra_runinfo.py
?? bin/validate_fastq_pair.awk
?? docs/refactor_review.md
?? lib/InputMetadata.groovy
?? lib/SraMetadata.groovy
?? modules/local/merge_sra_runs/
?? modules/local/sra_runinfo/
?? subworkflows/local/sra_acquisition/
```

## git diff --stat

Este comando muestra sólo cambios de archivos ya versionados; los archivos nuevos sin seguimiento se enumeran arriba y todavía no están en este stat. No se añadieron al índice para obtener la estadística.

```text
 .gitignore                                         |   10 +
 README.md                                          | 1187 ++------------
 assets/raw/test_metagenomic_samples.tsv            |    9 -
 bin/acquire_sra_sample.py                          |  700 ---------
 bin/check_samplesheet.py                           |  182 ---
 bin/manage_sra_checkpoints.py                      | 1619 --------------------
 bin/monitor_storage.py                             |   62 +-
 bin/resolve_sra_project.py                         | 1280 ----------------
 bin/spades_coassembly.py                           |   66 +-
 bin/summarize_resources.py                         |   45 +-
 conf/disk_efficient.config                         |   12 +-
 conf/modules.config                                |   21 +-
 conf/resources.config                              |   16 +-
 containers/sra-tools/Dockerfile                    |    3 +-
 containers/sra-tools/README.md                     |   17 +-
 containers/sra-tools/environment.yml               |    2 +-
 docs/architecture.md                               |  941 +-----------
 docs/compression_audit.tsv                         |   14 +-
 docs/original_resource_mapping.tsv                 |   10 +-
 metagenomics_pipeline.sh                           | 1156 +++-----------
 modules/core/sratools/environment.yml              |    2 +-
 modules/core/sratools/main.nf                      |   55 +-
 modules/core/sratools/meta.yml                     |   40 +-
 modules/local/check_samplesheet/environment.yml    |    5 -
 modules/local/check_samplesheet/main.nf            |   33 -
 modules/local/sra_preprocessing/environment.yml    |    5 -
 modules/local/sra_preprocessing/main.nf            |  146 --
 modules/local/sra_project/environment.yml          |    6 -
 modules/local/sra_project/main.nf                  |   92 --
 nextflow.config                                    |   13 +-
 nextflow_schema.json                               |  103 +-
 subworkflows/local/local_input/main.nf             |  115 +-
 .../local/sra_checkpoint_reconciliation/main.nf    |   22 -
 subworkflows/local/sra_global/main.nf              |  105 --
 subworkflows/local/sra_input/main.nf               |   73 +-
 subworkflows/local/sra_project_discovery/main.nf   |   51 -
 .../local/sra_sample_preprocessing/main.nf         |  112 --
 tests/fixtures/resource_accounting/human_trace.tsv |    3 -
 tests/fixtures/resource_accounting/raw_trace.tsv   |    3 -
 .../fixtures/resource_accounting/sra_manifest.tsv  |    4 -
 .../storage_usage_timeseries.tsv                   |    3 -
 .../resource_accounting/task_workdir_peaks.tsv     |    4 -
 .../task_workdir_timeseries.tsv                    |    5 -
 .../resource_accounting/trace_registry.tsv         |    4 -
 tests/fixtures/runinfo_invalid.csv                 |   10 -
 tests/fixtures/runinfo_selection.csv               |    9 -
 tests/fixtures/runinfo_valid.csv                   |    5 -
 tests/fixtures/sra_samples_invalid.tsv             |    9 -
 tests/fixtures/sra_samples_valid.tsv               |    2 -
 tests/test_check_samplesheet.py                    |  115 --
 tests/test_coverm_contig_awk.py                    |   95 --
 tests/test_coverm_mag_abundance_awk.py             |  104 --
 tests/test_genemarks2_rendering.py                 |  252 ---
 tests/test_launcher_lock_contract.py               |  118 --
 tests/test_monitor_storage.py                      |  250 ---
 tests/test_normalize_abundance.py                  |  268 ----
 tests/test_phylophlan_rendering.py                 |  270 ----
 tests/test_slurm_accounting.py                     |  102 --
 tests/test_sra_helpers.py                          | 1387 -----------------
 tests/test_static_contracts.py                     |  637 --------
 tests/test_summarize_resources.py                  |  264 ----
 workflows/metagenomics.nf                          |   51 +-
 62 files changed, 609 insertions(+), 11695 deletions(-)
```
