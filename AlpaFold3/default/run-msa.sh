#!/bin/sh
#PBS -l select=1:ncpus=64:mpiprocs=1:ompthreads=64
#PBS -l walltime=165:00:00

if [ ! -z "${PBS_O_WORKDIR}" ]; then
  cd "${PBS_O_WORKDIR}"
fi

module -s purge
. /apl/alphafold/af300-2_init.sh

RUNAF3=/apl/alphafold/run-af-300-2.sh

# notes:
# - YOU NEED TO OBTAIN MODEL PARAMETER FILE BY YOURSELF.
# - https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
# - actual options to run_alphafold.py can be displayed with "--dryrun -e".
# - input file for inference is generated in output_dir as
#   (output_dir)/(name)/(name)_data.json
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
   -j ./nqr.json \
   -o ./af_output \
   -M
