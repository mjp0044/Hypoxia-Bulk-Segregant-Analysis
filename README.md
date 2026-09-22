> **Note to self before pasting into GitHub:** the image paths below mirror this project's existing folder
> structure (e.g. `QTL figures and results/Gprime.jpg`). Upload those folders alongside `README.md` in the repo,
> or edit the paths to match wherever you land the figures. Also double-check the "Study design" paragraph below —
> I inferred the Pcrit/bulk-segregant framing from file and variable names; adjust the wording to match your
> actual methods description.

# Bulk-Segregant QTL Mapping of Hypoxia Tolerance in *Tigriopus californicus*

Bulk segregant analysis (BSA) / QTL-seq mapping of loci underlying hypoxia tolerance in the intertidal copepod
*Tigriopus californicus*, using pooled whole-genome sequencing of phenotypically extreme bulks from F30 hybrids
between two source populations (SD and SH). Full analysis code: [`QTLseq_BSA_script.R`](QTLseq_BSA_script.R).

## Study design

Two DNA pools ("Top" and "Bottom" bulks, sampled from opposite tails of the hypoxia-tolerance phenotype
distribution) were sequenced and their allele counts merged into a single table of biallelic SNPs. Genome-wide
divergence between bulks was quantified across 12 chromosomes using [QTLseqr](https://github.com/bmansfeld/QTLseqr)'s
G′ statistic (750 kb sliding windows), and significant regions were annotated with [SnpEff](https://pcingola.github.io/SnpEff/)
and cross-referenced against gene models, a prior RNA-seq differential-expression time series, and curated
functional gene sets.

## QTL mapping

G′ was calculated genome-wide and significant QTLs called at α = 0.01:

![G-prime genome-wide scan](QTL%20figures%20and%20results/Gprime.jpg)

| Chromosome | QTL | Region (Mb) | Length (Mb) | # SNPs | Max G′ | Mean q-value |
|---|---|---|---|---|---|---|
| Chr_11 | 1 | 11.86 – 14.95 | 3.10 | 15,510 | 4.52 | 0.0020 |
| Chr_12 | 2 | 0.02 – 1.98 | 1.96 | 6,601 | **5.62** | 0.0019 |
| Chr_12 | 3 | 6.96 – 8.23 | 1.27 | 3,560 | 4.17 | 0.0030 |
| Chr_3 | 4 | 4.95 – 5.38 | 0.43 | 2,220 | 3.57 | 0.0086 |
| Chr_3 | 5 | 6.07 – 7.21 | 1.14 | 4,718 | 3.92 | 0.0049 |
| Chr_8 | 6 | 9.92 – 11.14 | 1.22 | 6,534 | 3.85 | 0.0042 |
| Chr_9 | 7 | 13.08 – 14.45 | 1.37 | 4,505 | 3.77 | 0.0057 |

**7 significant QTLs** spanning **5 of 12 chromosomes**, mapping to **721 candidate genes**.

<details>
<summary>Supplementary QC/diagnostic plots</summary>

| SNP density across the genome | Chromosomes carrying significant peaks | ΔSNP-index with confidence intervals |
|---|---|---|
| ![nSNPs](QTL%20figures%20and%20results/nSNPs.jpg) | ![Gprime peaks](QTL%20figures%20and%20results/Gprime_chroms_with_peaks.jpg) | ![deltaSNP](QTL%20figures%20and%20results/deltaSNP.jpg) |

</details>

## Zooming in on a candidate gene

Each QTL was annotated down to the SNP level with SnpEff, letting individual candidates be inspected in genomic
context — here, a missense-carrying gene (*cyp6a14* / TCAL_08505) on Chr_11, with ΔSNP-index plotted against the
local gene model:

![cyp6a14 gene track](Gene%20Tracks/cyp6a14%20only%20gene%20track%20with%20delta%20snp.jpg)

(Several other candidates — *Tret1*, *Lcc2*, *GFPT1*, *Gnpda1*, *ATPsynCf6*, *Cha*, *Oxa1l*, *mesh* — were
profiled the same way; see the [`Gene Tracks/`](Gene%20Tracks) folder.)

## Modeling gene-level SNP architecture

Beyond calling QTLs, the 721 candidate genes were characterized with a series of generalized linear models, each
matched to the response variable's distribution:

- **Quasi-Poisson regression** (with `offset(log(gene_length_kb))`) for SNP counts falling in different genomic
  contexts (exons, introns, upstream regions, SnpEff impact categories) — the offset converts raw counts into an
  implicit per-kb rate while keeping the response on its natural count scale.
- **Gamma regression** (log link) for continuous, strictly-positive metrics like gene length and per-gene peak G′.
- **Linear models** on variance-stabilized (sqrt- or log-transformed) ratios where a Gaussian residual structure
  was more appropriate.

Each model compared two groupings of interest: whether a candidate gene overlapped an independent RNA-seq
differential-expression (DE) time series, and which source population's allele rose in frequency across a QTL
(SD- vs. SH-biased). Fitted group means and 95% CIs from these models are overlaid on raw per-gene distributions:

| DE-overlap comparison | Population-bias comparison |
|---|---|
| ![Ridgeline, DE overlap](Time%20Series%20overlap%20comparisons/ridgeline_column_plot.png) | ![Ridgeline, population bias](Pop%20bias%20comparisons/ridgeline_column_plot_pb.png) |

Three representative model fits, in full:

<details>
<summary><b>In-gene SNP density vs. RNA-seq DE overlap</b> — quasi-Poisson, gene-length offset</summary>

DE-overlapping genes carry ~21% more SNPs per gene than expected from gene length alone.

```r
summary(ingene.mod)
# n_in_gene ~ In_time_series + offset(log(gene_length_kb)), family = quasipoisson

Coefficients:
                   Estimate Std. Error t value Pr(>|t|)    
(Intercept)         1.61493    0.02403  67.212  < 2e-16 ***
In_time_seriesTRUE  0.18906    0.05517   3.427 0.000645 ***
---
(Dispersion parameter for quasipoisson family taken to be 10.68484)

    Null deviance: 7884.2  on 720  degrees of freedom
Residual deviance: 7763.5  on 719  degrees of freedom
```

</details>

<details>
<summary><b>Intronic SNP burden vs. population bias</b> — quasi-Poisson, gene-length offset</summary>

SH-biased QTL genes carry ~16% fewer intronic SNPs (length-adjusted) than SD-biased genes.

```r
summary(inintron.mod.pb)
# n_introns ~ bias + offset(log(gene_length_kb)), family = quasipoisson

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.33184    0.03087  43.148  < 2e-16 ***
biasSH_bias -0.16900    0.06110  -2.766  0.00582 ** 
---
(Dispersion parameter for quasipoisson family taken to be 11.2823)

    Null deviance: 8903.0  on 720  degrees of freedom
Residual deviance: 8814.3  on 719  degrees of freedom
```

</details>

<details>
<summary><b>Peak QTL signal (max G′) vs. population bias</b> — Gamma GLM, log link</summary>

SH-biased QTLs have ~9% lower peak G′ than SD-biased QTLs — the strongest population-bias effect in the dataset.

```r
summary(gprime.mod.pb)
# max_Gprime ~ bias, family = Gamma(link = "log")

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.395559   0.004461  312.82   <2e-16 ***
biasSH_bias -0.088774   0.008022  -11.07   <2e-16 ***
---
(Dispersion parameter for Gamma family taken to be 0.009911458)

    Null deviance: 7.7406  on 720  degrees of freedom
Residual deviance: 6.5406  on 719  degrees of freedom
AIC: 628.88
```

</details>

## Which population's allele rose in frequency?

Within each QTL, allele frequency trajectories were tracked separately for the SD- and SH-source alleles across
bulks, letting each QTL be classified as SD- or SH-biased:

![Allele frequency by QTL and bulk](Allele%20Frequency%20Plots/AF_by_QTL_and_bulk_plot.jpg)

Of the 721 candidate genes, **498 (69%) fall in SD-biased QTLs** and **223 (31%) fall in SH-biased QTLs**
(Chr_3/11/12 vs. Chr_8/9) — a split that turned out to be predictive of both QTL signal strength and intronic SNP
burden (models above).

## Do QTL genes line up with an independent RNA-seq study?

Candidate genes were cross-referenced against a prior RNA-seq differential-expression (DE) time series. **101 of
721 candidate genes (14%)** were also DE in that study — a modest, non-significant enrichment over the
genome-wide background rate (Fisher's exact test, OR = 1.18, p = 0.13). Genes that *did* overlap, however, carry
significantly more SNPs per gene after adjusting for gene length (model above), along with a higher burden of
intronic and "modifier"-impact SNPs specifically.

## What are the QTL genes doing functionally?

GO enrichment (`topGO`, Fisher's exact test, BP ontology) was run on the 721 candidate genes and visualized with
`rrvgo` after semantic-similarity clustering:

![GO term treemap](TopGO%20and%20rrvgo%20results/rrvgo%20treemap%20all.jpg)

Targeted overlap tests against curated functional gene sets (candidate vs. genome-wide background, Fisher's exact
test) did not show significant enrichment for any single category tested, though point estimates trended in
different directions:

| Gene set | Odds ratio | p-value |
|---|---|---|
| Mitochondrial-targeted proteins | 0.83 | 0.43 |
| Cuticle / chitin processing | 1.07 | 0.75 |
| Glycolysis & related processes | 0.53 | 0.23 |

## Methods at a glance

- **QTL calling:** [`QTLseqr`](https://github.com/bmansfeld/QTLseqr) G′ statistic, 750 kb sliding windows, α = 0.01
- **Variant annotation:** [SnpEff](https://pcingola.github.io/SnpEff/)
- **Gene-level modeling:** GLMs (quasi-Poisson with length offsets, Gamma with log link) and linear models on
  transformed ratios, comparing candidate-gene subgroups (RNA-seq DE overlap, SD-/SH-population bias)
- **GO enrichment:** `topGO` (weight01 algorithm, Fisher's exact test) + `rrvgo` (Wang semantic similarity,
  custom `org.Tcalifornicus.eg.db` OrgDb) for visualization
- **Visualization:** `ggplot2`, `ggridges`, `gggenes`, `patchwork`

Full pipeline, including data wrangling, filtering thresholds, and every model/plot summarized above, is in
[`QTLseq_BSA_script.R`](QTLseq_BSA_script.R).
