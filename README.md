# XCIE-EM

X-chromosome inactivation escape and phasing from single cell sequencing data using an expectation-maximization algorithm.

### Glossary

- Haplotype: The two X chromosomes present in female cells. When parent of origin information (which chrX was inherited from which parent) is available the haplotypes are called `mat` and `pat` for maternally and paternally inherited chromosomes. When parent of origin is not available the haplotypes are called `hap1` and `hap2`.
- Variant phase: The variant phase of heterozygous variants means which allele is present in which chromosomes. Here we always refer to the location of the reference allele, so for example phase `mat` means that the reference allele is in the maternal haplotype and alternate allele is in the paternal haplotype.
- Active chrX: Which of the two chrX is active within a cell. The chrX which is not inactivated.
- Cell activity: Refers to the active chrX within a cell. Analogous to variant phase (tells which chromosome the ref allele is in), cell activity tells which chromosome is active in a cell.

## Compiling

```
git clone https://github.com/maickrau/XCIE-EM.git
cd XCIE-EM
git submodule update --init --recursive
make all
```

## Running

First align the single cell reads to your reference and call variants. Make sure that the cell barcodes are in the BAM tags `CB:Z:`. Then run:

```
XCIE-EM --input-bam sorted_alignments.bam --input-vcf sample_variants.vcf --output-prefix EM_output --exclude-PAR --exclude-XIST-grch38 --annotation-gff3 gene_annotation_file.gff3
```

This will count ref/alt coverage per heterozygous variant per cell, and then estimate the variant phases and cell activities. The output is written to files which start with `EM_output`. See the section [output](#output) for explanations of all the output files.

Alternatively you can use [scReadCounts](https://horvathlab.github.io/NGS/SCReadCounts/) to count the expression per heterozygous SNP per cell and then run:

```
XCIE-EM --input-screadcounts screadcounts.output.tsv --output-prefix EM_output --exclude-PAR --exclude-XIST-grch38 --annotation-gff3 gene_annotation_file.gff3
```

The parameters `--exclude-PAR` `--exclude-XIST-grch38` are not mandatory but they improve the results by removing the PAR and XIST regions, which are difficult to phase from chrX inactivation. For other references you can also use `--exclude-XIST-grch37` or `--exclude-XIST-chm13` or manually with `--exclude-region`.

### Using trio data

If you have trio phasing data for some of the variants of the sample, you can use that to force the trio phased variants into correct phase. The variants which are not phased in the trio phased data will still be phased by the EM algorithm. Even a small amount of trio phased variants will help guide the phasing of the rest of the variants. The trio phasing data is inputed with the parameter `--force-phase trio_phase_file.tsv`. The required format of the trio phase file is one line per variant, with two columns for variant name and phase separated by tab, for example:

```
X:12345:A:G	mat
X:12356:T:C	pat
X:12456:T:G	pat
```

Variant name should be in format `X:(position):(ref_allele):(alt_allele)` and phase should be one of `mat` or `pat`, which means that the reference allele comes from maternal or paternal chromosome respectively. Alternately if the variants are phased without knowing their parent of origin, phase can be given as `hap1` or `hap2`, meaning the reference allele is in first or second haplotype respectively. The two formats should not be mixed so all phases should be either `mat`/`pat` or `hap1`/`hap2`. Variants outside of chrX and variants which don't match any variants in the input are allowed and silently ignored.

## Output

- `output.variants.tsv`: Description of variant phases, and phase confidence. The phase columns refers to the phase of the reference allele, so eg. phase `mat` means the reference allele comes from the maternal chromosome.
- `output.cells.tsv`: Description of cell activities, and activity confidence. The cell activity confidence is estimated only from variants which are estimated to be non-escape. The active column refers to the active chrX, so eg. active `mat` means that the maternally inherited chrX is active and paternally inherited chrX is inactive.
- `output.cells.withescapevariants.tsv`: Description of cell activities, and activity confidence. The cell activity confidence is estimated from all variants, which biases the cell activity confidence estimate downwards for cells which have high expression in escape variants. Not recommended to be used unless you know what you are doing.
- `output.pseudobulk.variants.confidence2.tsv`: The pseudobulk expression of each variant. The expression is split by haplotype and by active/inactive expression. Only variants and cells which are confidently phased are included.
- `output.pseudobulk.variants.confidence0.tsv`: The pseudobulk expression of each variant. The expression is split by haplotype and by active/inactive expression. All variants and cells are included, even those whose phasing is uncertain.
- `output.pseudobulk.genes.confidence2.tsv`: The pseudobulk expression of protein coding chrX genes. The expression is split by haplotype and by active/inactive expression. This sums over the variants in `output.pseudobulk.variants.confidence2.tsv`. Requires the parameter `--annotation-gff3`.
- `output.pseudobulk.genes.confidence0.tsv`: The pseudobulk expression of protein coding chrX genes. The expression is split by haplotype and by active/inactive expression. This sums over the variants in `output.pseudobulk.variants.confidence0.tsv`. Requires the parameter `--annotation-gff3`.
- `output.preprocessed_matches.tsv`: A file with the cell vs variant matches preprocessed into a suitable format.

### Phase confidence and activity confidence

The phase confidence of variants in the output files means how confident the probabilistic model is that the given phase is correct. It refers to the difference in log-likelihood between the current phase and the other phase, that is, `confidence = ln [ P(observed data | current variant phase) / P(observed data | flipped variant phase) ]`. The activity confidence of cells is calculated similarly from the log-likelihood difference.

### Haplotype names

If parent of origin is given from trio data then the variants and cells are phased using the terms `mat` and `pat`, where variant phase `mat` means the reference allele is maternal and cell activity `mat` means the maternal chrX is active. If no parent of origin is available then the haplotypes are called `hap1` and `hap2`, where again variant phase `hap1` means that the reference allele comes from the first haplotype and cell activity `hap1` means that the chrX represented by the first haplotype is active. That is, the terms `hap1` and `hap2` are consistent between the variant phase and cell activity.

### Pseudobulk p-values

The pseudobulk variant and gene files have columns `pvalue_Xiover10_unadjusted` and `pvalue_Xiunder10_unadjusted`. These are the p-values of binomial test for checking if the fraction of inactive expression is above / below 10%. The column `Xi` also shows the fraction of inactive expression. The p-values are **not adjusted for multiple testing**. You should correct the p-values for multiple testing.

## Parameters

- `--input-bam` and `--input-vcf`: Count SNPs from the given BAM file and use those as the input. Recommended. Requires that the cell barcode is in BAM tag `CB:Z:`.
- `--input-screadcounts`: Take the counts from scReadCounts as input. Recommended.
- `--input-preprocessed-table`: Take the counts from a preprocessed table.
- `-o, --output-prefix`: Output prefix.
- `--force-phase`: Use variant phasing from trio data to force the phase for some variants. This will also help guide the phasing of other variants.
- `--annotation-gff3`: A gff3 annotation file. Used for counting pseudobulk expression per gene.
- `--EM-noise-magnitude`: Noise added to the EM algorithm to avoid local minima. Default 20 usually works well.
- `--EM-noise-decay`: Decay speed of noise in EM algorithm. Default 0.95 usually works well.
- `--EM-random-seed`: Random seed used for EM initialization.
- `--EM-num-runs`: Number of initializations used for EM. The best result is picked out of all runs.
- `--exclude-regions`: Exclude some regions. The reason to exclude regions is that some regions (eg XIST) have higher expression from the inactive chrX, which violates the fundamental assumptions used in the EM algorithm.
- `--exclude-PAR`: Exclude the PAR region. Recommended since the PAR region usually has low phasing accuracy.
- `--exclude-XIST-grch38 --exclude-XIST-grch37 --exclude-XIST-chm13`: Exclude the XIST and TSIX genes. Recommended since XIST has higher expression from the inactive chrX. You should pick the correct one based on your reference.
