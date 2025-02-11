# **JuncSeq**

JuncSeq is a bioinformatics tool designed to detect and analyse aberrant splicing events using splice junction data. It enables researchers to:
1. Identify aberrant splicing in a predefined list of genes.
2. Detect splice junction abnormalities near specific genetic variants.

## Installation

You can install JuncSeek directly from GitHub:

```sh
pip install git+https://github.com/CCICB/juncseek.git
```

Alternatively, you can clone the repository and install it in editable mode:

```sh
git clone https://github.com/CCICB/juncseek.git
cd juncseek
pip install --editable .
```

## Usage

Once installed, you can run `juncseek` from the command line as follows:

```sh
juncseek [-h] --patient-list PATIENT_LIST --gtf GTF [--gene-list GENE_LIST] [--vcf VCF] --output OUTPUT
```

### **Arguments:**
- `--patient-list PATIENT_LIST` : Path to the patient list file (only needs paths as input).
- `--gtf GTF` : Path to the GTF annotation file.
- `--gene-list GENE_LIST` *(optional)* : Path to the gene list file.
- `--vcf VCF` *(optional)* : Path to the VCF file containing variant data.
- `--output OUTPUT` : Path to the output directory.
