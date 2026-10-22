# SVShare

A tool for comparing, filtering, and annotating structural variants across samples.

## Development Setup

A [Conda](https://docs.conda.io/) environment configuration is provided for convenience.

```bash
conda env create -f environment.yml
conda activate svshare
```

## Workflow

SVShare is organized into two main stages: **Analyze** and **Report**.

### Analyze

1. **Provide BAM or CRAM files** for each sample.
2. **Run Sniffles2 and cuteSV** on each sample to identify structural variants.
3. **Use Jasmine to group calls of the same SV** across both callers and all samples, and save the results in a SQLite database.
4. **Compare structural variants against gnomAD and CoLoRSdb** to determine how common or rare they are in the population.
5. **Use BEDTools with gene annotations** to identify affected genes and genomic regions.
6. **Save the complete analysis results.**

### Report

Filter the analysis results and write the SVs that pass as an HTML report or a TSV table. SVShare provides default filters to highlight rare structural variants shared across samples. Users can adjust these filters for their own analysis.

## Command-Line Interface

Run SVShare from the repo root as `python3 -m svshare <command> [options]`, where the command is `analyze` or `report`.

### Analyze

```bash
python3 -m svshare analyze \
  --samples sample1.bam sample2.bam sample3.bam \
  --reference GRCh38.fa \
  --output results/ \
  --gnomad gnomad.v4.1.sv.sites.vcf.gz \
  --colorsdb CoLoRSdb.GRCh38.v1.2.0.pbsv.jasmine.vcf.gz \
  --gencode gencode.v50.annotation.gtf.gz
```

All three are optional, but recommended.

- `--gnomad`: the gnomAD v4.1 SV sites VCF, for looking up how common each SV is.
  - Download it from the [gnomAD downloads page](https://gnomad.broadinstitute.org/downloads), with its `.tbi` index in the same directory.
- `--colorsdb`: the CoLoRSdb GRCh38 `pbsv.jasmine` VCF, also for looking up how common each SV is.
  - Download it from [CoLoRSdb v1.2.0 on Zenodo](https://zenodo.org/records/14814308), with its `.tbi` index in the same directory.
  - It's built from long reads, so it includes SVs that gnomAD's short reads miss, such as many insertions and SVs in repetitive regions.
- `--gencode`: GENCODE's human gene annotation, for finding the genes each SV affects.
  - Download `gencode.v50.annotation.gtf.gz` from the [GENCODE human releases page](https://www.gencodegenes.org/human/).
  - Each gene found is also given a readable name from [HGNC](https://www.genenames.org/), a summary of what it does from [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/), and dosage sensitivity scores from the [ClinGen Dosage Sensitivity Map](https://dosage.clinicalgenome.org/). These come from `svshare/data/gene_info.tsv`, which is included with SVShare; to refresh it, run `python3 scripts/build_gene_info.py`.

### Report

Running `report` without additional options writes an HTML report of the SVs that pass the default filters below:

```bash
python3 -m svshare report \
  --results results/
```

Choose the output with `--format`:

- `html` (default): written to `<results>/report.html`.
- `tsv`: a table with one row per SV, written to `<results>/report.tsv`.

Pass `--output` to write either one elsewhere. Both include only the SVs that pass the filters.

Available filter options:

```text
--caller-support {both,sniffles2,cutesv,any}   Default: any
--sv-type {DEL,DUP,INS,INV,BND,any}            Default: any
--max-gnomad-frequency <number>                Default: 0.01
--min-samples <number>                         Default: 2
--genes <gene1> <gene2> ...
--gene-file <file>
--region <genomic-region>
```

For longer gene lists, users can provide a text file:

```bash
python3 -m svshare report \
  --results results/ \
  --gene-file genes.txt
```

The gene list file should contain **one gene symbol per line** with no header:

```text
BRCA1
BRCA2
TP53
PALB2
```
