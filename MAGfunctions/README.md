This folder contains data and analyses to produce a functional summary of 40 MAGs. we usued [gapseq](https://github.com/jotech/gapseq) and [antiSMASH](https://antismash.secondarymetabolites.org/#!/download) to predict complete functions and biosynthetic gene clusters (BGCs), respectively. Results from gapseq were curated either with the help of pathway description file from [(Weiss et al 2021)](https://www.nature.com/articles/s41396-021-01153-z) or manually by using [MetaCyc database](https://metacyc.org/META/organism-summary) .

To produce the combined data for the final figure, run scripts in following order:

1.antismash.R

2.gapseq.R

3.gapseq_manual.R

4.combined_antismash_gapseq.R
