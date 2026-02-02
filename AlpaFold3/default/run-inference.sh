#!/bin/sh
#PBS -l select=1:ncpus=16:mpiprocs=1:ompthreads=16:ngpus=1
#PBS -l walltime=72:00:00

if [ ! -z "${PBS_O_WORKDIR}" ]; then
  cd "${PBS_O_WORKDIR}"
fi

module -s purge
. /apl/alphafold/af300-2_init.sh

RUNAF3=/apl/alphafold/run-af-300-2.sh

# notes:
# - YOU NEED TO OBTAIN MODEL PARAMETER FILE BY YOURSELF.
# - input json file for inference is the output of data-pipeline (MSA) run
# - actual options to run_alphafold.py can be displayed with "--dryrun -e".
# - GPU is required for inference.
#
# options:
#   -j <input json>, --json_path=<input json>
#     Path to input json file (single input file).
#   -i <input directory>, --input_dir=<input directory>
#     Path to directory which contains input json files (multiple input files).
#   -o <output directory>, --output_dir=<output directory>
#     Path to output directory.
#   -m <model directory>, --model_dir=<model directory>
#     Path to model for inference (THIS IS NOT PROVIDED BY RCCS).
#   -M, --msa, --msaonly
#     Do only data pipeline (MSA).
#   -I, --inference, --inferenceonly
#     Skip MSA and do only inference. (GPU required)
#   -e, --echo
#     Echo actual options to run_alphafold.py
#   --dryrun
#     Do not run actually.
#
#   Other run_alphafold.py options may also be accepted.

$RUNAF3 \
   -j ./af_output/nqr_nadh_uq9/nqr_nadh_uq9_data.json \
   -o ./af_output \
   -m ~/AlphaFold3/models \
   -I
