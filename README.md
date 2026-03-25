# LPSS_BLAST
Code for the paper: Structure-based proteome mining reveals latent symbiotic nitrogen-fixation modules in rice.

## Installation
### Getting started
```bash
git clone https://github.com/wangge-1996/LPSS_BLAST.git
cd LPSS_BLAST
```
### Install Dependencies
```yaml
dependencies:
  - python=3.9
  - pandas>=2.2.3
  - biopython>=1.84
  - pip>=24.3.1
```
## Usage
Prior to running the pipeline, you must configure Tm-vec according to its official documentation, and save the downloaded pretrained model to the `model` folder. Official repo: https://github.com/tymor22/tm-vec
```bash
python LPSS_BLAST.py --query_gene query.fasta --database_gene target.fasta
```

## Citation
