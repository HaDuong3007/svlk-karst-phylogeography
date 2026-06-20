# svlk-karst-phylogeography

Data, scripts, and manuscript for:

> Ha Duong *et al.* (submitted). Nuclear–chloroplast marker discordance and cross-border geographic structuring in two karst-endemic plant genera (*Begonia* sect. *Coelocentrum*, Begoniaceae; *Hemiboea*, Gesneriaceae) of the Sino-Vietnamese Limestone Belt: an exploratory analysis. *Journal of Asia-Pacific Biodiversity*.

---

## Overview

The Sino-Vietnamese Limestone Belt (SVLB) is a near-continuous karst corridor straddling the Vietnam–China border (Ha Giang, Cao Bang, Lang Son ↔ Guangxi, Yunnan, Guizhou). This repository contains everything needed to reproduce the maximum-likelihood phylogenetic analyses, AU topology tests, and figures reported in the manuscript.

Two plant genera were examined:
- ***Begonia* sect. *Coelocentrum*** (Begoniaceae) — ITS and rpl16, 23–24 taxa
- ***Hemiboea*** (Gesneriaceae) — ITS (32 taxa) and trnL-F (9 taxa)

All sequences are publicly available on GenBank. No new sequences were generated. All ML analyses used **IQ-TREE 3.1.1**.

---

## Repository layout

```
svlk-karst-phylogeography/
├── alignments/              MAFFT v7.526 alignments + specimen metadata
│   ├── begonia_svlb_ITS_aligned.fasta       (23 taxa, 788 bp)
│   ├── begonia_svlb_rpl16_aligned.fasta     (24 taxa, 1071 bp)
│   ├── begonia_svlb_metadata.csv
│   ├── hemiboea_main_ITS_aligned.fasta      (32 taxa, 895 bp)
│   ├── hemiboea_main_metadata.csv
│   ├── hemiboea_trnlf_aligned.fasta         (9 taxa, 893 bp)
│   └── hemiboea_trnlf_metadata.csv
│
├── trees/                   IQ-TREE 3.1.1 ML treefiles + full logs
│   ├── begonia_ITS_fullsupport.treefile     (TN+F+R2, logL=−4064.06)
│   ├── begonia_ITS_fullsupport.iqtree
│   ├── begonia_rpl16_fullsupport.treefile   (K3Pu+F+G4, logL=−2504.96)
│   ├── begonia_rpl16_fullsupport.iqtree
│   ├── hemiboea_ITS_triple.treefile         (HKY+F+G4, logL=−3523.21)
│   ├── hemiboea_ITS_triple.iqtree
│   ├── hemiboea_trnlf_triple.treefile       (F81+F, logL=−1431.02)
│   └── hemiboea_trnlf_triple.iqtree
│   (Node labels: SH-aLRT/aBayes/UFBoot)
│
├── au_test/                 AU topology test for Begonia ITS vs rpl16
│   ├── begonia_ITS_shared.fasta             (23-taxon shared ITS alignment)
│   ├── begonia_rpl16_shared.fasta           (23-taxon shared rpl16 alignment)
│   ├── both_topologies.nwk                  (ITS and rpl16 ML topologies)
│   ├── au_ITS.iqtree                        (rpl16 topology vs ITS data)
│   └── au_rpl16.iqtree                      (ITS topology vs rpl16 data)
│
├── figures/                 Manuscript figures (PDF + PNG, 600 dpi)
│   ├── fig1_begonia_trees.pdf/png    Figure 1  — Begonia ITS (A) + rpl16 (B)
│   ├── fig2_begonia_map.pdf/png      Figure 2  — Begonia distribution map
│   ├── fig3_hemi_trnlf.pdf/png       Figure 3  — Hemiboea trnL-F tree (9 taxa)
│   └── figS1_hemi_ITS.pdf/png        Figure S1 — Hemiboea ITS tree (32 taxa)
│
├── scripts/                 Full analysis pipeline (01–32)
│   ├── 01_genbank_inventory.py       GenBank inventory retrieval
│   ├── 02–20 ...                     Data curation, alignment, QC
│   ├── 21_begonia_viz.R              Begonia trees + map (intermediate)
│   ├── 27_fetch_trnlf.py             Fetch trnL-F sequences
│   ├── 28_build_concat.py            Concatenation setup
│   ├── 30_hemi_multimarker_viz.R     Hemiboea multi-marker figures (intermediate)
│   ├── 31_final_figures.R     ← Main script: all 4 manuscript figures
│   └── 32_au_test.py          ← AU test setup and execution
│
├── manuscript/
│   ├── manuscript_draft_v4_JAPB.md          (submitted text)
│   └── supplementary_table_S1.md            (all 42 GenBank accessions)
│
├── README.md                (this file)
├── .zenodo.json             (Zenodo deposit metadata)
├── .gitignore
└── LICENSE                  (CC-BY 4.0)
```

> **Note on script numbering:** Scripts follow the pipeline execution order. Numbers 22–26 and 29 correspond to intermediate curation and QC steps handled interactively and are not represented as standalone scripts. Script 31 (`31_final_figures.R`) produces all four manuscript figures and is the primary visualisation script for end users.

---

## Reproducing the analyses

### Prerequisites

| Tool | Version used | Install |
|------|-------------|---------|
| IQ-TREE | 3.1.1 | `conda install -c bioconda iqtree` |
| MAFFT | 7.526 | `conda install -c bioconda mafft` |
| R | ≥ 4.1 | https://www.r-project.org |
| Python | ≥ 3.9 | https://www.python.org |

R packages: `ggplot2`, `ggtree`, `treeio`, `ape`, `dplyr`, `patchwork`, `sf`, `rnaturalearth`

```r
install.packages(c("ggplot2", "ape", "dplyr", "patchwork"))
BiocManager::install(c("ggtree", "treeio"))
install.packages(c("sf", "rnaturalearth", "rnaturalearthdata"))
```

### ML trees (IQ-TREE 3.1.1)

Run from the repository root. Pre-computed treefiles are in `trees/`; these commands reproduce them from scratch.

```bash
# Begonia ITS  (TN+F+R2)
iqtree3 -s alignments/begonia_svlb_ITS_aligned.fasta \
        -m TN+F+R2 -B 1000 --bnni --alrt 1000 --abayes \
        --prefix trees/begonia_ITS_fullsupport

# Begonia rpl16  (K3Pu+F+G4)
iqtree3 -s alignments/begonia_svlb_rpl16_aligned.fasta \
        -m K3Pu+F+G4 -B 1000 --bnni --alrt 1000 --abayes \
        --prefix trees/begonia_rpl16_fullsupport

# Hemiboea ITS  (HKY+F+G4)
iqtree3 -s alignments/hemiboea_main_ITS_aligned.fasta \
        -m HKY+F+G4 -B 1000 --bnni --alrt 1000 --abayes \
        --prefix trees/hemiboea_ITS_triple

# Hemiboea trnL-F  (model selected by ModelFinder BIC → F81+F)
iqtree3 -s alignments/hemiboea_trnlf_aligned.fasta \
        -m MFP -B 1000 --bnni --alrt 1000 --abayes \
        --prefix trees/hemiboea_trnlf_triple
```

### AU topology test

```bash
# Requires IQ-TREE, R/ape, and Python ≥ 3.9
python3 scripts/32_au_test.py
# Outputs written to au_test/
```

The script auto-detects `iqtree3`/`iqtree2` and `Rscript` from PATH. Pre-computed AU test logs are already in `au_test/`.

### Figures

```bash
# All four manuscript figures (Fig. 1, Fig. 2, Fig. 3, Fig. S1) → figures/
Rscript scripts/31_final_figures.R
```

The script auto-resolves paths for either the full pipeline directory structure or this clean repository layout, and outputs both PDF (vector) and PNG (600 dpi) for each figure.

---

## Key results

| Dataset | Model | logL | Key node | Support (SH-aLRT/aBayes/UFBoot) |
|---------|-------|------|----------|----------------------------------|
| *Begonia* ITS | TN+F+R2 | −4064.06 | {*B. cavaleriei* + *B. longa* + *B. mashanica*} | 70.3/0.43/57 |
| *Begonia* rpl16 | K3Pu+F+G4 | −2504.96 | {*B. cavaleriei* + *B. masoniana*} | 79.4/0.992/73 |
| *Hemiboea* ITS | HKY+F+G4 | −3523.21 | {*H. longgangensis* + *H. mollifolia* + *H. longzhouensis*} | 85.0/0.984/76 |
| *Hemiboea* trnL-F | F81+F | −1431.02 | {*H. follicularis* + *H. longgangensis* + *H. longzhouensis* + *H. rubribracteata*} | 83.4/0.995/73 |

**AU test (Begonia ITS vs rpl16):**

| Test | ΔlogL | p-AU |
|------|-------|------|
| rpl16 topology on ITS data | 149.96 | 2.13 × 10⁻⁶³ — **rejected** |
| ITS topology on rpl16 data | 273.32 | 1.78 × 10⁻⁵ — **rejected** |

---

## Sequence provenance

All sequences from GenBank. Full provenance (accession, depositor, original publication, locality, habitat) in `manuscript/supplementary_table_S1.md`.

| Source paper | Accessions | Marker |
|-------------|-----------|--------|
| Chung *et al.* 2014 *Bot. Stud.* 55:1 | KF636xxx, KF707xxx | *Begonia* ITS, rpl16 |
| Möller *et al.* 2009 *Am. J. Bot.* 96:2108 | FJ501356–357 | *Hemiboea* ITS |
| Möller *et al.* 2011 *Plant Syst. Evol.* 292:223 | HQ632xxx | *Hemiboea* ITS, trnL-F |
| Li XQ *et al.* 2022 *Proc. R. Soc. B* 289:20211308 | OK322552–553 | *Hemiboea* ITS |
| Liang *et al.* 2025 *Ecol. Evol.* 15:e72330 | PV439816–PV448327 | *Hemiboea* ITS |
| Huang J *et al.* 2017 (direct submission) | KY607413 | *Hemiboea* trnL-F |

---

## Citation

If you use these data or scripts, please cite the manuscript (above) and the original sequence sources listed in `manuscript/supplementary_table_S1.md`.

**Software citations:**
- IQ-TREE 3: Wong *et al.* 2025 (submitted); Minh *et al.* 2020 *Mol. Biol. Evol.* 37:1530
- MAFFT: Katoh & Standley 2013 *Mol. Biol. Evol.* 30:772
- ModelFinder: Kalyaanamoorthy *et al.* 2017 *Nat. Methods* 14:587
- UFBoot2: Hoang *et al.* 2018 *Mol. Biol. Evol.* 35:518

---

## License

This repository is released under the [Creative Commons Attribution 4.0 International (CC-BY 4.0)](LICENSE) licence. You are free to share and adapt the material for any purpose, provided appropriate credit is given.

---

## Contact

Ha Duong
