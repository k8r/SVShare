# SVShare

A tool for comparing, filtering, and annotating structural variants across samples.

## Workflow

SVShare is organized into three main stages: **Analyze**, **Filter**, and **Report**.

### Analyze

1. **Provide BAM or CRAM files** for each sample.
2. **Run Sniffles2 and cuteSV** on each sample to identify structural variants.
3. **Use Truvari to identify variants supported by both callers.**
4. **Use Jasmine to identify shared structural variants across samples.**
5. **Compare shared structural variants against gnomAD** to determine how common or rare they are in the population.
6. **Use BEDTools with gene annotations** to identify affected genes and genomic regions.
7. **Save the complete analysis results.**

### Filter

Use an existing SVShare analysis to apply user-selected filters, such as caller support, SV type, gnomAD frequency, number of samples, genes, or genomic region.

### Report

Generate a formatted report from either the complete or filtered SVShare results.

## Command-Line Interface

### Analyze

```bash
svshare analyze \
  --samples sample1.bam sample2.bam sample3.bam \
  --reference GRCh38.fa \
  --output results/
```

### Filter

```bash
svshare filter \
  --results results/ \
  --caller-support both \
  --max-gnomad-frequency 0.01 \
  --sv-type DEL \
  --min-samples 2 \
  --genes BRCA1 BRCA2 TP53
```

Possible filter options include:

```text
--caller-support {both,sniffles2,cutesv,any}
--sv-type {DEL,DUP,INS,INV,BND,any}
--max-gnomad-frequency <number>
--min-samples <number>
--genes <gene1> <gene2> ...
--gene-file <file>
--region <genomic-region>
```

For longer gene lists, users can provide a text file:

```bash
svshare filter \
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

### Report

```bash
svshare report \
  --results results/ \
  --output report.html
```
