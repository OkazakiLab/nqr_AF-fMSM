This repository contains **scripts and associated data for constructing and analyzing Markov state models (MSMs) from molecular dynamics (MD) trajectories**, as used in the analysis of the following preprint:

- Preprint: https://doi.org/10.64898/2026.01.28.702413

## Directory structure
```text
nqr_AF-fMSM
├── AlphaFold3/                # job scripts & input files for AlphaFold3
│   ├── default/
│   │   ├── nqr.json
│   │   ├── run_msa.sh         # sample script for MSA preparation
│   │   └── run_inference.sh   # sample script for structural inference
│   ├── af3_mmm/               # sample file for MSA subsampling
│   │                            (cf. https://github.com/OkazakiLab/af3_mmm)
│   └── template_usage/        # sample file for template usage
├── MSM_FESred/                # MD data for MSM under reduced 2Fe-2S conditions
│   ├── msm_data/
│   ├── AF3/                   # MD data initiated from AF3 structures
│   ├── round0/                # MD data initiated from experimental or AF2 structures
│   ├── round1/                # MD data in adaptive round 1
│   ├── round2/                # MD data in adaptive round 2
│   └── round3/                # MD data in adaptive round 3
├── MSM_FESox/                 # MD data for MSM under oxidized 2Fe-2S conditions
├── scripts/
│   ├── mdanalysis.py          # functions for MD data analysis
│   ├── msm.py                 # functions for MSM analysis
│   ├── 0_mdanalysis_exe.py    # executable script for MD data analysis
│   ├── 1_tica_exe.py          # executable script for data preparation for MSM analysis
│   ├── 2_msm_exe.py           # executable script for MSM analysis
│   ├── input_features.csv
│   └── exp_data/
└── README.md
```

## Computational environments
### MD trajectory analysis
- Python 3.13.2  
- MDAnalysis v2.9.0  
- NumPy v2.2.4  
- etc.

### Markov state modeling
- Python 3.9.18  
- PyEMMA v2.5.12  
- NumPy v2.0.2  
- etc.

> Note: Separate Python environments are used for MD analysis and MSM construction.

## How to run the analysis
The analysis should be performed in the following order:

0. Feature extraction from MD trajectories
python3 0_mdanalysis_exe.py

1. tICA and data preparation for MSM analysis
python3 1_tica_exe.py

2. MSM analysis and visualization
python3 2_msm_exe.py

