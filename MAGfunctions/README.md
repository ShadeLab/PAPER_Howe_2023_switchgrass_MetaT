* This folder contains data and analyses to produce a functional summary of 40 MAGs. we usued [gapseq](https://github.com/jotech/gapseq) and [antiSMASH](https://antismash.secondarymetabolites.org/#!/download) to predict complete functions and biosynthetic gene clusters (BGCs), respectively. Results from gapseq were curated either with the help of pathway description file from [(Weiss et al 2021)](https://www.nature.com/articles/s41396-021-01153-z) or manually by using [MetaCyc database](https://metacyc.org/META/organism-summary) .

To produce the combined data for the final figure, run scripts in following order:

1.antismash.R

2.gapseq.R

3.gapseq_manual.R

4.combined_antismash_gapseq.R

* Associated methods in text:

Biosynthetic gene clusters (BGC) were predicted by antiSMASH (v6.0) (ref) and further annotated by Big-SCAPE (v1.1.0). 

We used gapseq (ref) to predict completed metabolic pathways in our focal MAGs. We used the ‘find -p all –b 200’ option to search for pathways against the MetaCyc database. We filtered out incomplete pathways and remaining pathways were grouped into broader categories using MetaCyc classification and manually curated to highlight only relevant. These groups were defined by potential involvement: i) plant (using plant metabolites/cell components), ii) phytohormone (known/potential involvement in phytohormone homeostasis), iii) stress (e.g., drought, reactive oxygen species), and iv) general (pathways that utilize potential plant derived products). Furthermore, we also manually searched for genes/pathways that were missed by gapseq but are known to be related to adaptation to plant-associated lifestyle, including secretion systems, oxidation of trace gases (H2 and CO), oxalate degradation, and phytohormone production/degradation.  
