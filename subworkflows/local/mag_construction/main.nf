#!/usr/bin/env nextflow

include { ASSEMBLY } from './assembly/main'
include { BINNING } from './binning/main'
include { MAG_REFINEMENT } from './mag_refinement/main'

// TODO Phase 2+: orchestrate ASSEMBLY, BINNING, and MAG_REFINEMENT.
// This skeleton intentionally produces no biological output.
workflow MAG_CONSTRUCTION {}
