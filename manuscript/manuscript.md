# Decoding post-transcriptional regulatory networks by RNA-linked CRISPR screening in human cells

Patrick J. Nugent<sup>1,2</sup>, Heungwon Park<sup>1</sup>, Cynthia L. Wladyka<sup>3</sup>, Jamie Yelland<sup>1</sup>, Sayantani Sinha<sup>4</sup>, Katharine Y. Chen<sup>1,2</sup>, Christine Bynum<sup>1,5</sup>, Grace Quarterman<sup>1,5</sup>, Stanley C. Lee<sup>4</sup>, Andrew C. Hsieh<sup>3,6</sup>, Arvind Rasi Subramaniam<sup>1,7,†</sup>

<sup>1</sup> Basic Sciences Division and Computational Biology Section of the Public Health Sciences Division, Fred Hutchinson Cancer Center, Seattle WA, USA  
<sup>2</sup> Molecular and Cellular Biology Graduate Program, University of Washington, Seattle WA, USA  
<sup>3</sup> Human Biology Division, Fred Hutchinson Cancer Center, Seattle WA, USA  
<sup>4</sup> Translational Science and Therapeutics Division, Fred Hutchinson Cancer Center, Seattle WA, USA  
<sup>5</sup> Department of Biology, Spelman College, Atlanta GA, USA  
<sup>6</sup> Department of Medicine and Department of Genome Sciences, University of Washington, Seattle WA, USA  
<sup>7</sup> Department of Biochemistry and Department of Genome Sciences, University of Washington, Seattle WA, USA

<sup>†</sup> Corresponding author: A.R.S: <rasi@fredhutch.org>

## Abstract

RNAs undergo a complex choreography of metabolic processes that are regulated by thousands of RNA-associated proteins. Here we introduce ReLiC, a scalable and high-throughput RNA-linked CRISPR approach to measure the responses of diverse RNA metabolic processes to knockout of 2,092 human genes encoding all known RNA-associated proteins. ReLiC relies on an iterative strategy to integrate genes encoding Cas9, sgRNAs, and barcoded reporter libraries into a defined genomic locus. Combining ReLiC with polysome fractionation reveals key regulators of ribosome occupancy, uncovering links between translation and proteostasis. Isoform-specific ReLiC captures differential regulation of intron retention and exon skipping by SF3b complex subunits. Chemogenomic ReLiC screens decipher translational regulators upstream of mRNA decay and identify a role for the ribosome collision sensor GCN1 during treatment with the anti-leukemic drug homoharringtonine. Our work demonstrates ReLiC as a powerful framework for discovering and dissecting post-transcriptional regulatory networks in human cells.

## Introduction

After transcription, RNAs undergo several metabolic events such as splicing, editing, localization, translation, and decay inside cells. These RNA processes are executed by ribonucleoprotein complexes composed of RNA-binding proteins (RBPs), adapter proteins, and regulatory factors. Over 2,000 human genes encode proteins that are part of ribonucleoprotein complexes<sup>[1](#ref-Gerstberger2014),[2](#ref-Hentze2018)</sup>. Mutations in RNA-associated proteins occur in many human diseases such as cancer, neurodegeneration, and developmental disorders<sup>[3](#ref-Gebauer2021)</sup>.

Crosslinking-based biochemical approaches can identify RBP-RNA interactions<sup>[2](#ref-Hentze2018)</sup> but not their functional effect on RNA metabolism. RBPs can increase, decrease or leave unchanged metabolic events on their target RNA depending on their affinity, location, and other associated factors<sup>[4](#ref-VanNostrand2020),[5](#ref-Ray2013)</sup>. Many RBPs also associate with multiple protein complexes and participate in several distinct RNA metabolic events<sup>[6](#ref-Schneider-lunitz2021)</sup>. Conversely, protein factors that do not directly bind RNA can still affect RNA metabolism by regulating the interactions between RNAs and RBPs or by controlling the cellular level and activity of RBPs<sup>[7](#ref-England2022)</sup>.

Unbiased genetic screening can identify cellular factors regulating RNA metabolism, but are limited in their current form. CRISPR screens using indirect phenotypes such as cell growth and fluorescent protein levels are difficult to engineer and interpret for many RNA metabolic events<sup>[8](#ref-Przybyla2022),[9](#ref-Genolet2022)</sup> due to potential false positives<sup>[10](#ref-Gonatopoulos-Pournatzis2018),[11](#ref-Scarborough2021)</sup> and genetic compensatory mechanisms<sup>[12](#ref-El-Brolosy2019)</sup>. CRISPR perturbations followed by pooled single cell RNA sequencing can capture RNA phenotypes such as steady-state levels, polyadenylation status, and decay rates<sup>[13](#ref-Replogle2022)–[15](#ref-Xu2024)</sup>. But these transcriptome-wide approaches have limited flexibility to study different types of RNA processes, show bias towards highly expressed RNAs, and remain costly and labor intensive to scale beyond a few dozen perturbations. Thus, genome scale screening approaches to identify the RNA-centric functions of human proteins, that also have the flexibility to capture diverse RNA metabolic events, are highly desirable.

## Results

### Development of RNA-linked CRISPR screening in human cells

We reasoned that combining CRISPR-based perturbations with barcoded RNA readouts provides a general approach to study the genetic control of RNA processes in human cells. Indeed, RNA interference screens in human cells<sup>[16](#ref-Herholt2018)</sup> and CRISPR interference screens in *S. cerevisiae*<sup>[17](#ref-Muller2020),[18](#ref-Alford2021)</sup> have linked perturbations to barcoded transcriptional readouts. Lentiviral delivery, commonly used for CRISPR screening in human cells, scrambles sgRNA-barcode linkages due to template switching during reverse transcription<sup>[19](#ref-Sack2016)</sup> and results in variable expression of RNA barcodes due to random genomic integration<sup>[20](#ref-Ellis2005)</sup>. To avoid these limitations, we used an iterative, Bxb1-mediated site-specific integration strategy to stably express doxycycline-inducible SpCas9 (Cas9 hereafter), sgRNAs, and barcoded RNA reporters from a defined genomic locus ([Fig. 1a](figures/fig1.png)). Using an *EYFP* reporter, we confirmed its uniform expression and Cas9-mediated depletion in HEK293T ([Fig. 1b](figures/fig1.png)) and U2OS cells ([Extended Data Fig. 1a,b](figures/figs1.png)).

To identify post-transcriptional regulators, we targeted 2,092 human genes encoding known RNA-associated proteins<sup>[1](#ref-Gerstberger2014),[2](#ref-Hentze2018)</sup> ([Fig. 1c](figures/fig1.png)). We used a dual sgRNA design with random N<sub>20</sub> barcodes in a modular vector allowing insertion of arbitrary RNA reporters ([Extended Data Fig. 1c](figures/figs1.png)). Our final library targeted 2,190 genes with 4 sgRNA pairs per gene, including non-targeting and essential gene-targeting controls (Supplementary Table 1). We linked the N<sub>20</sub> barcodes to sgRNAs by paired-end deep sequencing of the cloned plasmid library. We then integrated this library into HEK293T cells, and counted barcodes in the genomic DNA and transcribed RNA by deep sequencing ([Fig. 1c](figures/fig1.png), [Extended Data Fig. 1d](figures/figs1.png)). We recovered a median of 8 barcodes per sgRNA pair (henceforth referred to as sgRNAs) with at least one barcode for 99% of sgRNAs and 100% of all genes ([Extended Data Fig. 1e](figures/figs1.png)), thus capturing the diversity of our input library.

To test whether sgRNA-linked barcodes capture fitness effects, we counted barcodes in genomic DNA and mRNA after Cas9 induction (Supplementary Table 6). Barcode counts showed little systematic change 5 days after Cas9 induction ([Fig. 1d](figures/fig1.png), left panels). However, on days 13 and 21 after Cas9 induction, barcode counts for a subset of sgRNAs were strongly depleted in both the genomic DNA and mRNA in a highly correlated manner ([Fig. 1d](figures/fig1.png), middle and right panels). The magnitude of depletion was correlated across distinct barcode sets for each gene ([Extended Data Fig. 1f](figures/figs1.png)), indicating the assay’s technical reproducibility. Barcodes in both genomic DNA and mRNA corresponding to annotated essential genes (n = 745) were depleted at days 13 and 21 relative to other genes targeted by our library (n = 1401, [Fig. 1e](figures/fig1.png)). Thus, RNA-linked CRISPR captures both the identity and fitness effect of genetic perturbations solely from sequencing of barcodes in mRNA and genomic DNA.

### ReLiC identifies regulators of mRNA translation

We first applied ReLiC to study translation, an RNA process that is not directly accessible in existing CRISPR screening methods. We combined ReLiC with polysome profiling<sup>[21](#ref-Warner1963)–[23](#ref-Gierer1963)</sup> to separate mRNAs based on their ribosome occupancy. We used a *β*-globin reporter<sup>[24](#ref-Zhang1998)</sup> as a model of a well-translated mRNA ([Fig. 2a](figures/fig2.png)), inserted random barcodes into its 3<sup>′</sup> UTR, and confirmed that over 75% of the *β*-globin mRNA was in polysome fractions ([Fig. 2b](figures/fig2.png)). We then cloned the *β*-globin reporter into our ReLiC plasmid library, integrated the plasmids into HEK293T cells, induced Cas9 for 7 days, and fractionated cell lysates ([Fig. 2a](figures/fig2.png)). After counting sgRNA-linked barcodes in pooled fractions, we used MAGeCK<sup>[25](#ref-Li2014)</sup> to identify sgRNAs that significantly altered the ratio of barcode counts between heavy (H) or light (L) polysomes and monosomes (M) (Supplementary Table 7). Polysome to monosome ratios for individual sgRNAs were highly reproducible (r=0.92 and 0.80 for H/M and L/M, respectively) between replicate experiments ([Extended Data Fig. 2a](figures/figs2.png)). We used an FDR threshold of 0.05 and a minimum of three concordant sgRNAs for calling gene hits that altered polysome to monosome ratios (Supplementary Table 8).

304 and 207 gene hits decreased heavy polysome to monosome and light polysome to monosome ratios, respectively ([Fig. 2c](figures/fig2.png)). 37 gene hits increased heavy polysome to monosome ratio, while 2 increased light polysome to monosome ratio ([Fig. 2c](figures/fig2.png)). 176 gene hits overlapped between the heavy and light polysome ratios, indicating a largely concordant effect of gene knockouts on different polysome fractions. The skewed distribution of gene hits with more perturbations decreasing ribosome occupancy likely arises from the efficient translation of *β*-globin mRNA in unperturbed cells ([Fig. 2b](figures/fig2.png)). Consistent with heavy polysome fractions containing better-translated mRNAs, heavy polysome to monosome ratios were more sensitive to perturbations with more gene hits and larger effect sizes than light polysome to monosome ratios ([Fig. 2c](figures/fig2.png)). We therefore focused on heavy polysome to monosome ratios for further analyses.

Gene hits that decreased polysome to monosome ratios were highly enriched for cytoplasmic ribosomal proteins and ribosome biogenesis factors ([Extended Data Fig. 2b](figures/figs2.png)). On average, knockout of large ribosomal proteins and biogenesis factors decreased polysome to monosome ratios more than knockout of small ribosomal proteins and biogenesis factors ([Fig. 2d](figures/fig2.png)). Sequencing of barcodes in supernatant fractions revealed that knockout of small ribosomal proteins and biogenesis factors led to increased barcode representation in the supernatant relative to monosomes, while knockout of large ribosomal proteins and biogenesis factors had the opposite effect ([Fig. 2e](figures/fig2.png)). This observation is consistent with small ribosomal subunit depletion preventing association between mRNAs and ribosomes, while large ribosomal subunit depletion still allowing scanning by one or more small ribosomal subunits. Indeed, barcode ratios between supernatant and polysome revealed comparable effects of small and large ribosomal perturbations ([Extended Data Fig. 2c](figures/figs2.png)).

Depletion of most translation initiation factors also decreased heavy polysome to monosome ratios, but their effects were generally smaller than the effect of ribosomal protein depletion ([Fig. 2f](figures/fig2.png)). Some initiation factor subunits not classified as hits (open circles, [Fig. 2f](figures/fig2.png)) still had multiple sgRNAs that decreased heavy polysome to monosome ratio, but either fell just below our gene-level FDR threshold (EIF4G1), or did not meet our stringent criterion of 3 distinct sgRNAs with significant effects (EIF2S2, EIF4E). In the case of the 12-subunit EIF3 and associated EIF3J, the seven subunits A,B,C,D,E,G, and I called as hits are the same ones that severely reduce polysome to monosome ratio and fitness when depleted by siRNA in HeLa cells<sup>[26](#ref-Wagner2016)</sup>. While most canonical initiation factors emerged as hits upon comparison of polysome fraction to either the monosome or the supernatant fractions, a few genes such as XRN1, DDX6, and SNTB1 emerged as hits only in the polysome to supernatant comparison ([Extended Data Fig. 2d](figures/figs2.png)). Aminoacyl-tRNA synthetase knockouts had mild and variable effects on ribosomal occupancy ([Fig. 2f](figures/fig2.png)), presumably reflecting a balance between their direct effect on elongation and indirect effect on initiation through GCN2 and EIF2*α* phosphorylation<sup>[27](#ref-Darnell2018)</sup>.

We identified several knockouts outside the core translation machinery with decreased polysome to monosome ratio ([Fig. 2f](figures/fig2.png)). Subunits of the CCR4-NOT complex (CNOT1, CNOT2, CNOT3, and CNOT7), which has been implicated in a wide range of RNA metabolic processes<sup>[28](#ref-Collart2016)</sup>, emerged as hits in our screen. Knockout of proteasome and TRiC chaperonin subunits led to substantially reduced polysome to monosome ratios, comparable in magnitude to knockout of core translation initiation factors ([Fig. 2f](figures/fig2.png)). We validated these effects by creating individual knockout lines, and observed reduced bulk polysome to monosome ratios after 3 days of proteasomal or chaperonin depletion ([Extended Data Fig. 2e](figures/figs2.png)). These complexes did not arise as hits only due to their essentiality since knockout of other essential cellular complexes such as RNA polymerase II and splicing factor 3a/b (SF3) did not reduce polysome to monosome ratio ([Fig. 2f](figures/fig2.png)). While neither the proteasome nor the TRiC chaperonin complex has been directly associated with translational regulation, they play a critical role in maintaining cellular proteostasis by coordinating their activities with translation<sup>[29](#ref-Harper2016),[30](#ref-Frydman2012)</sup>. Our results suggest a reciprocal regulation of translation in response to changes in proteasomal and chaperonin capacity.

### Relation between ribosome occupancy and growth fitness

Ribosome occupancy on mRNAs is often correlated with cellular growth rate, with slower growth accompanied by lower polysome to monosome ratio across different growth conditions and organisms<sup>[26](#ref-Wagner2016),[31](#ref-Balakrishnan2022),[32](#ref-Metzl-Raz2017)</sup>. Our measurements of both depletion of barcodes and their polysome distribution across thousands of gene perturbations allowed testing the generality of the relationship between ribosome occupancy and growth. Across all perturbations, decrease in polysome to monosome ratio was positively correlated with barcode depletion in both mRNA and genomic DNA but had a wide distribution ([Extended Data Fig. 2f](figures/figs2.png)). However, gene groups corresponding to different molecular complexes had characteristically distinct relationships between ribosome occupancy and barcode depletion. Perturbing ribosomal proteins and biogenesis factors resulted in the largest decrease in polysome to monosome ratio relative to fitness, followed by EIF3 ([Fig. 2g](figures/fig2.png), [Extended Data Fig. 2g](figures/figs2.png)). Perturbing proteasomal subunits produced a smaller but still significant decrease in ribosome occupancy, while perturbing RNA polymerase II subunits did not alter ribosome occupancy despite their significant effects on barcode depletion ([Fig. 2g](figures/fig2.png), [Extended Data Fig. 2g](figures/figs2.png)). Hence, the coupling between growth rate and ribosome occupancy in human cells is not invariant across all perturbations, but depends on the pathway or the molecular process that leads to growth limitation.

### Regulatory factors that increase ribosome occupancy

We next examined the small group of gene knockouts that increased the heavy polysome to monosome ratio ([Fig. 2h](figures/fig2.png), brown triangles). Knockout of EEF2, EIF5A, and EEF1A increased polysome to monosome ratio, consistent with their role in promoting translation elongation. Intriguingly, the ribosome-associated quality control factor ASCC3 was the top gene hit for increased heavy polysome to monosome ratio (log<sub>2</sub>H/M = 0.62, FDR = 1e-4). Since ASCC3 is involved in splitting stalled ribosomes on mRNAs<sup>[33](#ref-Juszkiewicz2020)</sup>, its presence here suggests that even well-translated mRNAs such as *β*-globin undergo some degree of ribosome stalling and quality control. Further, knockout of the ribosome collision sensor ZNF598, which acts upstream of ASCC3<sup>[33](#ref-Juszkiewicz2020)</sup>, also increased ribosome occupancy (log<sub>2</sub>H/M = 0.19, FDR = 0.06, p = 0.007). Knockout of METAP2, which removes methionine from the N-terminus of nascent polypeptides, increased ribosome occupancy (log<sub>2</sub>H/M = 0.22, FDR = 0.001, p = 3e-4), pointing to an effect of nascent peptide processing on the kinetics of mRNA translation.

Finally, we asked whether differential effects of gene perturbations on ribosome occupancy as measured by polysome to monosome ratios are reflected in their cellular transcriptional response. Using a genome-scale Perturb-seq dataset<sup>[13](#ref-Replogle2022)</sup>, we correlated and clustered the transcriptional profiles of translation factor perturbations that had concordant or discordant effects on ribosome occupancy ([Fig. 2i](figures/fig2.png)). Perturbations with concordant effects on ribosome occupancy ([Fig. 2h](figures/fig2.png)) did not show a higher correlation with each other than with perturbations with discordant effects on ribosome occupancy. For example, depletion of METAP2 and EIF2S1 (EIF2*α*), which are direct interactors<sup>[34](#ref-Datta1988)</sup>, had a markedly higher correlation in their transcriptional responses even though these gene knockouts had discordant effects on ribosome occupancy ([Fig. 2h](figures/fig2.png)). Thus, the effects of gene perturbations on ribosome occupancy measured by ReLiC are distinct from their downstream transcriptional responses.

### Isoform-specific splicing screens using ReLiC

Existing screening approaches to study RNA splicing require careful design of fluorescent protein reporters<sup>[10](#ref-Gonatopoulos-Pournatzis2018)</sup> and can result in high false positive and negative rates<sup>[11](#ref-Scarborough2021)</sup>. We reasoned that ReLiC will allow us to directly measure the ratio of different splice isoforms carrying the same barcode, thereby capturing the effect of the linked sgRNA perturbation on splicing. To test this idea, we used the same *β*-globin reporter as in our translation screen ([Fig. 3a](figures/fig3.png)), and confirmed the canonical isoform with three exons as the most abundant one ([Fig. 3b](figures/fig3.png)). We then performed three isoform-specific screens for regulators that increase intron 1 retention (i12), intron 2 retention (i23), or exon 2 skipping (e13) ([Fig. 3a](figures/fig3.png)). After harvesting RNA 1, 3, 5 and 7 days after Cas9 induction, we selectively amplified each isoform along with the barcode ([Fig. 3c](figures/fig3.png)). We measured the ratio of barcode counts between each isoform and the total RNA pool and used an FDR threshold of 0.05 and a minimum of three concordant sgRNAs for calling gene hits (Supplementary Table 8). For all three isoforms, the number of gene hits progressively increased with longer duration of Cas9 induction ([Extended Data Fig. 3a](figures/figs3.png)). Fewer gene knockouts increased the exon 2-skipped isoform in comparison to the two intron-retained isoforms at all time points ([Extended Data Fig. 3a](figures/figs3.png)). Effect sizes of gene hits were reproducible across distinct barcode sets for each gene ([Extended Data Fig. 3b](figures/figs3.png)) and specific to each isoform ([Extended Data Fig. 3c](figures/figs3.png)).

The three isoform-specific screens identified both common and unique sets of gene hits that were evident by automated gene ontology analysis ([Fig. 3d](figures/fig3.png)) and by manual inspection ([Fig. 3e](figures/fig3.png)). Gene hits in the two intron retention screens were dominated by core spliceosome components and splicing-associated factors (yellow circles and triangles, [Fig. 3e](figures/fig3.png)). Spliceosome hits were distributed throughout the splicing cycle starting from the tri-snRNP complex that is required to form the catalytically active spliceosome and included members from most known spliceosomal subcomplexes<sup>[35](#ref-Wahl2015)</sup>. Our screen also identified *trans* regulators of spliceosomal function such as CDK11B – a recently identified activator of the SF3b complex<sup>[36](#ref-Hluchy2022)</sup>, and BRF2 – an RNA polymerase III subunit required for transcription of U6 snRNA.

Retention of intron 1 was promoted by an additional group of gene knockouts that were enriched for mRNA translation and nuclear RNA exosome factors (red and brown triangles, [Fig. 3e](figures/fig3.png)). Loss of ribosomal proteins and translation factors might inhibit nonsense-mediated decay (NMD) of the intron 1-retained isoform, while the effect of nuclear RNA exosome components might be indirect through their role in ribosome biogenesis. While retention of either intron 1 or intron 2 will generate a premature termination codon (PTC), only the intron 1-retained isoform will have a splice junction and an associated exon-junction complex (EJC) downstream of the PTC, which is a well-known trigger for NMD<sup>[37](#ref-Cheng1994)</sup>. Consistent with a role for NMD, EJC components (MAGOH, EIF4A3, RBM8A) and RNA export factors (NCBP1, NCBP2) emerged as hits only in the intron 1 retention screen ([Fig. 3e](figures/fig3.png)), and several NMD factors such as UPF2 and SMG1 increased intron 1 retention even though they fell below our FDR threshold for calling hits (Supplementary Table 8).

### Differential effects of SF3b complex subunits on splicing

In contrast to intron retention, perturbations increasing exon 2 skipping were enriched for a narrow set of splicing factors. Components of the U2 snRNP, most notably several members of the SF3 complex, were among the top hits (purple squares, [Fig. 3e](figures/fig3.png)), suggesting that their depletion allows some degree of splicing but impairs the correct selection of splice sites. This is consistent with the subtle alterations in exon skipping caused by disease-causing mutations in the SF3b complex<sup>[38](#ref-Darman2015)</sup>. Exon 2 skipping was also promoted by perturbing components involved in nuclear protein import (green squares, [Fig. 3e](figures/fig3.png)), presumably through their effect on nuclear import of U2 snRNP proteins after their synthesis in the cytoplasm. Perturbing individual components of the 7-subunit SF3b complex<sup>[39](#ref-Cretu2016)</sup> had distinct effects on exon skipping and intron retention ([Fig. 3f](figures/fig3.png)), even though all 7 subunits are essential for cell growth ([Extended Data Fig. 3d](figures/figs3.png)). Exon 2 skipping was greatly increased upon loss of the subunits SF3B1, SF3B2, SF3B3, SF3B5, slightly increased by loss of SF3B7, and unaffected by loss of SF3B4 and SF3B6 ([Fig. 3f](figures/fig3.png)). Intron 2 retention was increased by loss of SF3B6 and SF3B7, while intron 1 retention was increased by loss of SF3B1, SF3B2, and SF3B5 ([Fig. 3f](figures/fig3.png)). By contrast, loss of the activating helicase AQR increased the retention of both introns 1 and 2 (brown markers, [Fig. 3f](figures/fig3.png)).

We next examined how the differential effects of SF3b subunit depletion on *β*-globin reporter splicing extend to endogenous mRNAs. To this end, we generated HEK293T cell lines with the subunits SF3B5 and SF3B6, which affected distinct splicing events in our screen, individually depleted. We also targeted AQR, a top hit in both our intron retention screens, as a positive control and included a non-targeting control sgRNA (FLUC). We performed RNA-seq 4 days after Cas9 induction to identify endogenous splicing events that are particularly sensitive to the respective genetic perturbations. Loss of SF3B5 increased skipping of 45 annotated cassette exons by 10% or higher ([Fig. 3g](figures/fig3.png)). Loss of SF3B6 or AQR affected the skipping of less than 10 cassette exons at the same effect size, while all three splicing factors increased aberrant retention of a similar number of distinct introns ([Fig. 3g](figures/fig3.png)). Interestingly, for genes such as RPL24 and RPL41, increased intron retention and exon skipping upon SF3B5 loss occurred at distinct splice sites within the same transcriptional unit ([Fig. 3h,i](figures/fig3.png)). In summary, the differential effects of SF3b subunits on splicing of the *β*-globin reporter extend to endogenous mRNAs with a subset of SF3b subunits playing a more prominent role in regulating exon skipping.

### ReLiC screen for regulators of mRNA quality control

We reasoned that sequencing mRNA barcodes using ReLiC provides a general approach to identify regulators of mRNA quality control pathways independent of their effect on protein expression. To test this idea, we modified our *β*-globin reporter to add a premature termination codon (PTC) at position 39 in the second exon ([Fig. 4a](figures/fig4.png)), which is known to trigger NMD<sup>[24](#ref-Zhang1998)</sup>. At steady state, mRNA levels of the PTC-containing reporter were strongly reduced relative to a reporter with a normal termination codon (NTC, [Extended Data Fig. 4a](figures/figs4.png)). To measure mRNA effects specific to the PTC and NTC reporters, we combined our ReLiC-RBP library with a dual barcoding strategy<sup>[17](#ref-Muller2020)</sup> to normalize barcode counts for the reporter of interest relative to that of the mCherry-puro selection marker within each cell ([Fig. 4a](figures/fig4.png)). We harvested RNA 7 days after Cas9 induction and counted mRNA barcodes for the PTC and NTC *β*-globin reporters and the mCherry-puro marker.

We identified 90 gene hits (FDR \< 0.05, 3 sgRNAs with concordant effects) whose knockout increased levels of the PTC reporter relative to the mCherry-puro marker ([Fig. 4b](figures/fig4.png), [Extended Data Fig. 4c](figures/figs4.png)). We did not get any hits that increased mRNA levels of the NTC reporter, as expected from its higher mRNA stability ([Fig. 4b](figures/fig4.png), [Extended Data Fig. 4c](figures/figs4.png)). Several core components of the NMD pathway (UPF1, UPF2, SMG1, SMG5, SMG7, ETF1) were among the gene hits for the PTC reporter, indicating our ability to identify NMD-specific factors ([Fig. 4b](figures/fig4.png), pink circles). Other NMD-associated factors such as SMG6 and EIF4A3 fell just below the FDR threshold but still significantly (MAGeCK P-value \< 0.05) increased mRNA levels of the PTC reporter. Remarkably, a large proportion of gene hits for the PTC reporter encoded core factors involved in various steps of mRNA translation ([Fig. 4b](figures/fig4.png), squares, triangles, and diamonds; [Extended Data Fig. 4b](figures/figs4.png)). These included both small and large ribosomal proteins, ribosome biogenesis factors, translation initiation factors, and aminoacyl-tRNA synthetases. EIF2, EIF2B, and EIF3 subunits but not EIF4F subunits emerged as hits ([Extended Data Fig. 4d](figures/figs4.png)), despite having similar effects on growth fitness ([Extended Data Fig. 4e](figures/figs4.png)). Further, acute chemical inhibition of EIF4E and EIF4A showed only a modest effect on PTC reporter levels in comparison with SMG1 inhibition ([Extended Data Fig. 4a](figures/figs4.png)). Thus, our results suggest a limited *in vivo* role for EIF4F compared to EIF2, EIF2B, and EIF3 in regulating NMD of our *β*-globin reporter.

### Chemical and genetic modifier screens using ReLiC

Our NMD screen also identified gene hits involved in ER and mitochondrial homeostasis ([Fig. 4b](figures/fig4.png), x markers). This is consistent with phosphorylation of EIF2*α* upon perturbing ER and mitochondrial homeostasis leading to NMD inhibition<sup>[40](#ref-Wang2011)</sup>. To identify regulators of NMD acting through EIF2*α* phosphorylation, we adapted ReLiC to perform a chemical modifier screen using the small molecule ISRIB that renders translation insensitive to EIF2*α* phosphorylation<sup>[41](#ref-Sidrauski2013)</sup>. After inducing Cas9 for 6 days, we treated a ReLiC cell pool expressing the PTC reporter with ISRIB or DMSO for 48 hours, harvested RNA, and counted barcodes. We identified 30 gene knockouts (FDR \< 0.01) that decreased mRNA levels of the PTC reporter upon ISRIB treatment relative to the DMSO control ([Fig. 4c](figures/fig4.png)). These gene hits included several ER- and mitochondrially-localized proteins (Fig. 4c, x markers), consistent with their knockout inhibiting NMD through EIF2*α* phosphorylation.

Knockout of several aminoacyl-tRNA synthetases also decreased PTC reporter levels upon ISRIB treatment (Fig. 4c, diamonds), suggesting that their depletion inhibits NMD through phosphorylation of EIF2*α* rather than by decreasing translation elongation. To test this hypothesis, we performed a genetic modifier screen using ReLiC to deplete the EIF2*α* kinase GCN2, which is activated by uncharged tRNAs accumulating upon inhibition of aminoacyl-tRNA synthetases<sup>[42](#ref-Wek1995),[43](#ref-Dong2000)</sup>. We transduced the ReLiC cell pool with lentiviruses expressing either GCN2-targeting or non-targeting sgRNAs, induced Cas9 for 7 days, harvested RNA, and counted barcodes. Out of the 12 gene hits (FDR \< 0.01) with lower PTC reporter levels upon GCN2 depletion ([Fig. 4d](figures/fig4.png)), 10 were aminoacyl-tRNA synthetases (Fig. 4d, diamonds), thus confirming their action through GCN2-mediated EIF2*α* phosphorylation. Together, the above chemical and genetic modifier screens demonstrate the usefulness of ReLiC for dissecting the molecular pathways through which specific gene products regulate RNA processes.

### GCN1 regulates cellular responses to an anti-leukemic drug

Homoharringtonine (HHT) is an FDA-approved chemotherapeutic that targets the ribosome and is used to treat chronic myeloid leukemia<sup>[44](#ref-Gandhi2014)</sup>. HHT binds to the large ribosomal subunit to arrest initiating ribosomes at start codons and inhibit protein synthesis<sup>[45](#ref-Fresno1977),[46](#ref-Ingolia2011)</sup>, but how cells respond to this translational arrest is not well understood. Given ReLiC’s ability to identify regulators downstream of both mRNA translation and chemical perturbations, we sought to use this approach to probe the cellular response to HHT treatment. To this end, we performed ReLiC screens using an EYFP reporter ([Fig. 5a](figures/fig5.png)). After inducing Cas9 for 7 days, we treated the cell pool with 1 *μ*M HHT or DMSO for 6 hours before harvesting RNA and counting barcodes.

Unlike our previous ReLiC screens where we uncovered multiple gene hits and RNA metabolic pathways, a single gene, *GCN1*, emerged as a clear hit (FDR \< 0.05) whose knockout increased EYFP mRNA levels during HHT treatment ([Fig. 5b](figures/fig5.png)). GCN1 activates the kinase GCN2 to trigger EIF2*α* phosphorylation in response to amino acid limitation<sup>[47](#ref-Marton1993)</sup>. While GCN2 did not come up as a hit in our initial HHT screen, pharmacologic inhibition of GCN2 ablates the difference in EYFP mRNA levels between the GCN1 KO and control cells ([Extended Data Fig. 5a](figures/figs5.png)). GCN1 also binds collided ribosomes on mRNAs<sup>[48](#ref-Pochopien2021)</sup>, which can trigger both degradation of the nascent peptide and the mRNA<sup>[49](#ref-Oltion2023),[50](#ref-Muller2023)</sup>. Since ribosome collisions also trigger the ribotoxic stress response (RSR) through the kinase ZAK*α* that was not included in our original screen<sup>[51](#ref-Wu2020)</sup>, we measured p38 phosphorylation in wild-type and GCN1-depleted cells. HHT treatment increased p38 phosphorylation in GCN1-depleted cells ([Fig. 5d](figures/fig5.png)) in a ZAK*α*-dependent manner ([Extended data Fig. 5b](figures/figs5.png)), while wild-type cells did not show a corresponding increase. By contrast, treatment with the elongation inhibitor anisomycin triggered p38 phosphorylation in both wild-type and GCN1-depleted cells ([Fig. 5d](figures/fig5.png)).

Ribosome collisions induced by elongation inhibitors trigger upregulation of immediate early gene (IEGs) mRNAs<sup>[52](#ref-Sinha2020)</sup>. To test if GCN1 regulates IEGs during HHT treatment, we performed RNA-seq on wild-type and GCN1-depleted cells after HHT treatment ([Fig. 5e](figures/fig5.png)). HHT treatment for 6 hours caused widespread changes in mRNA levels in both wild-type and GCN1-depleted cells with \~225 up-regulated genes and \~450 down-regulated genes (\> 2-fold change, p \< 0.05). However, a small group of 60 genes, which included several immediate early genes such as *FOS*, *JUN*, and *MYC*, were differentially upregulated in GCN1-depleted cells relative to wild-type cells ([Fig. 5e](figures/fig5.png)). GCN1-depleted cells from the myeloid lineage AML cell line U937 also showed pronounced upregulation of IEGs upon HHT treatment ([Extended Data Fig. 5c](figures/figs5.png)). Further, pharmacologic inhibition of ZAK or GCN2 prevented IEG upregulation in GCN1-depleted cells during HHT treatment ([Extended Data Fig. 5d](figures/figs5.png)). The ZAK-dependent increased p38 signaling and IEG upregulation in GCN1-depleted cells suggest a role for GCN1 in mitigating ribosome collisions during HHT treatment.

To test whether ribosome collisions occur on endogenous mRNAs during HHT treatment, we first performed polysome fractionation from both wild-type and GCN1-depleted cells ([Extended Data Fig. 5e](figures/figs5.png)). Polysomes collapsed into monosomes after 1 hour of HHT treatment, and nuclease-resistant peaks, indicative of collided ribosomes, were of comparable intensity in both wild-type and GCN1-depleted cells. Additionally, ribosome profiling under the same conditions did not show obvious difference in average ribosome occupancy on mRNAs between wild-type and GCN1-depleted cells ([Fig. 5f](figures/fig5.png)). Thus, ribosome collisions do not occur during HHT treatment at a level that is detectable by bulk biochemical fractionation and do not alter global ribosome occupancy on mRNAs. Nevertheless, highly expressed immediate early genes such as *JUN* and *MYC* exhibit extensive ribosome density throughout their 5<sup>′</sup> UTR during HHT treatment ([Fig. 5g](figures/fig5.png), [Extended Data Fig. 5f](figures/figs5.png)). Furthermore, even in the absence of HHT, ribosomes initiate at multiple in-frame start codons on mRNAs of several immediate early genes such as *JUN*, *MYC*, and *JUND*<sup>[53](#ref-Hann1988)–[55](#ref-Gonzalez-Sanchez2024)</sup>. Together, these observations suggest that GCN1 senses collisions on these mRNAs between upstream scanning or elongating ribosomes and HHT-arrested initiating 80S ribosomes at downstream start codons.

## Discussion

Our study demonstrates ReLiC, an RNA-linked CRISPR approach for genetic dissection of diverse post-transcriptional processes in human cells. ReLiC enables measuring the effect of thousands of gene perturbations on mRNA translation, splicing, and decay – RNA processes that are not readily accessible to existing CRISPR screening methodologies. Our work reveals networks of molecular pathways, protein complexes, and individual proteins regulating these processes. These measurements are consistent with known molecular mechanisms and illuminate the complex interplay between different post-transcriptional regulatory events.

ReLiC reveals the role of essential pathways and genes in RNA metabolism even when their knockout is deleterious to cell growth. Chemical perturbations that abrogate protein expression can still be probed for their genetic dependencies, as seen from our identification of GCN1’s role during HHT treatment. ReLiC captures differential effects of perturbations within the same protein complex such as between members of the SF3b complex and between large and small ribosomal proteins. Unlike biochemical strategies, ReLiC identifies both direct effectors and indirect regulators, as exemplified by the identification of translation-related pathways across our screens for ribosomal occupancy, splicing, and mRNA decay.

ReLiC has complementary strengths and limitations for studying RNA metabolism compared to single cell approaches. ReLiC can be readily combined with bulk RNA sequencing-based readouts, thus rendering diverse RNA phenotypes such as localization, condensation, and editing amenable to genetic screening. Further, ReLiC can be used to probe rare events such as aberrant splicing that are difficult to capture using existing single cell approaches. Since ReLiC relies on targeted barcode sequencing, it can be executed at a genome scale in a cost-effective manner. Nevertheless, ReLiC depends on heterologously expressed reporters, and thus requires follow up studies on endogenous RNAs to establish transcriptome-wide significance. ReLiC also depends on non-lentiviral genomic integration to preserve linkage between sgRNA and reporter barcodes, which can limit its use to amenable cell types.

We anticipate that ReLiC can be extended to a broad range of biological settings, genetic perturbations, and RNA classes. Applying ReLiC to diverse cell types, cell states, and disease models will reveal differences in RNA metabolism that underlie cellular heterogeneity and disease progression. While we have used SpCas9 to induce gene knockouts, other effectors like base editors and prime editors can be readily incorporated into our modular workflow to identify the role of specific protein domains or regulatory elements on RNA metabolism with nucleotide resolution. Using non-coding, viral, and synthetic RNAs instead of mRNA reporters in ReLiC has the potential to unlock novel RNA regulatory mechanisms and therapeutic strategies. Finally, expanding ReLiC from our RNA interactome-focused library to all protein coding genes in the human genome will illuminate new interactions between RNA metabolism and other cellular processes.

## Acknowledgements

We thank members of the Subramaniam lab, the Basic Sciences Division, and the Computational Biology Program at Fred Hutch for assistance with the project and discussions. We thank Adam Geballe, Chris Lapointe, Akhila Rajan, and Brian Zid for feedback on the manuscript. This research was funded by NIH R35 GM119835 (A.R.S.), NSF MCB 1846521 (A.R.S.), NIH T32 GM008268 (P.J.N.), NIH R37 CA230617 (A.C.H.), NIH R01 CA276308 (A.C.H.), and NIH GM135362 (A.C.H.). This research was supported by the Genomics and Flow Cytometry Shared Resources of the Fred Hutch/University of Washington Cancer Consortium (P30 CA015704) and Fred Hutch Scientific Computing (NIH grants S10-OD-020069 and S10-OD-028685). The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.

## Author Contributions

P.J.N. designed research, performed experiments, analyzed data, and wrote the manuscript. H.P. performed experiments. C.L.W., J.Y., and A.C.H. assisted with polysome fractionation experiments. S.S and S.C.L. performed experiments on the U937 cell line. C.B., G.Q., and K.Y.C. performed gene ontology analyses. A.R.S. conceived the project, designed research, analyzed data, wrote the manuscript, supervised the project, and acquired funding.

## Competing Interests

None

## Captions

### Figure 1

**Development of RNA-linked CRISPR (ReLiC) screening.** **a.** *Strategy for genomic integration of Bxb1 landing pad, SpCas9, dual sgRNAs, and barcoded RNA reporters.* *attP* and *attP\** refer to orthogonal recombination sites for the Bxb1 integrase that differ by a single nucleotide mismatch and undergo recombination only with their corresponding *attB* and *attB\** sites. **b.** *Validation of Cas9 activity.* sgEYFP refers to an *EYFP*-targeting-sgRNA and sgCTRL is a non-targeting control. Histograms represent fluorescence of 10,000 cells as measured by flow cytometry. ‘Days post Cas9’ refers to days after addition of doxycycline to induce Cas9 expression. **c.** *Strategy for ReLiC sgRNA library design and validation.* **d.** *Correlated change in barcode frequency between genomic DNA and mRNA after Cas9 induction*. Each point corresponds to a gene knockout. Fold-changes are median-centered across all sgRNA pairs in the library. Gene level fold-changes are median values across all detected sgRNAs for each gene. *r* refers to Pearson correlation coefficient. **e.** *Essential gene knockouts are depleted in genomic DNA and mRNA after Cas9 induction.* Pan-essential genes from the DepMap database (n = 745) and remaining genes (n = 1401) are shown as separate histograms.

### Figure 2

**Polysome ReLiC identifies regulators of mRNA translation.** **a.** *Strategy for combining ReLiC and polysome fractionation.* **b.** *Reporter distribution across polysome fractions in unperturbed cells.* Points correspond to relative mRNA level in each fraction for distinct 3<sup>′</sup> UTR barcodes (n = 6) for the *β*-globin reporter. **c.** *Change in polysome to monosome ratio after 7 days of Cas9 induction.* Each point corresponds to a gene knockout. Horizontal axis indicates median of polysome to monosome ratios of barcode counts across all detected sgRNAs for each gene. Number of genes with FDR \< 0.05 and decreased or increased polysome to monosome ratio are indicated with *N* Individual gene hits are highlighted in dark grey triangles. All other genes are shown as light grey circles. **d.** *Change in polysome to monosome ratio for ribosomal protein and ribosome biogenesis genes.* Closed circles correspond to gene hits (FDR \< 0.05 with 3 or more concordant sgRNAs). **e.** *Change in monosome to supernatant ratio for ribosomal protein and ribosome biogenesis genes.* Left: equivalent to Fig. 2c but for monosome/supernatant. Ribosomal protein and ribosome biogenesis hits are highlighted. Right: equivalent to Fig. 2d but for monosome/supernatant. Vertical axes in c and e indicate P-values from a permutation test as calculated by MAGeCK. **f.** *Change in polysome to monosome ratio for protein groups and complexes.* Closed and open circles denote gene hits and non-hits similar to d. **g.** *Comparison of ribosome occupancy and mRNA depletion.* Points correspond to genes belonging to one of the highlighted groups. Shaded areas correspond to 95% confidence intervals for a linear fit. **h.** *Barcode ratios between polysome fractions for individual translation factors.* Each point corresponds to a distinct sgRNA and grey bars denote the median value across detected sgRNAs for that gene. **i.** *Correlation of expression profiles as measured by Perturb-seq<sup>[13](#ref-Replogle2022)</sup>.* *r* refers to Pearson correlation coefficient. EEF1A1, EEF1A2, and ZNF598 depletions did not induce strong transcriptional responses, so they are excluded from the visualization.

### Figure 3

**Isoform-specific splicing screens using ReLiC.** **a.** *Schematic of ReLiC splicing screens.* Location of RT primer and PCR primers used for PCR amplification of barcodes in each isoform are shown as black arrows. **b.** *Relative abundance of *β*-globin reporter splice isoforms as measured by RNA-seq.* Top panel shows RNA-seq coverage and bottom panel shows read counts mapping to each splice junction and intron. **c.** *Selective amplification of barcodes linked to splice isoforms.* Agarose gel lanes show RT-PCR products of expected size for the different isoforms. **d.** *Gene ontology analysis*. Fold enrichment of selected cellular processes and components 7 days after Cas9 induction. Dashes indicate GO terms with FDR \> 0.05. **e.** *Identity of splicing regulators.* Each point corresponds to a gene knockout. Isoform ratios are median values across all detected sgRNAs for each gene after median-centering across all sgRNAs in the library. Individual panels correspond to days after Cas9 induction (horizontal) and isoform screens (vertical). Marker shape denotes isoform identity and marker color denotes one of five highlighted gene groups. Genes with FDR \< 0.05 and belonging to one the highlighted groups are listed in the legend. **f.** *Relative reporter isoform levels upon SF3b complex perturbations.* AQR is shown as a positive control hit for intron retention. FDR \< 0.05 is indicated by large marker, and FDR ≥ 0.05 is indicated by small marker. **g.** *Change in endogenous splicing isoforms upon SF3b complex perturbations.* RNA-seq was performed 4 days after inducing Cas9. Change in intron retention or cassette exon skipping were calculated across all ENSEMBL-annotated transcripts, and ranked by decreasing magnitude of change with respect to the FLUC control sample. **h.** *Examples of endogenous isoform changes.* Changes in retained introns and skipped exons are highlighted in green and blue rectangles, respectively. Schematics at the bottom correspond to ENSEMBL isoforms with the highlighted retained intron and skipped exon events. **i.** *Quantification of isoform fraction* for the endogenous intron retention and exon skipping events in h.

### Figure 4

**Dissecting cotranslational quality control using chemogenomic ReLiC screening.** **a.** *Dual barcode strategy for measuring reporter mRNA levels.* Red octagons denote stop codons. **b.** *Gene hits from dual barcode NMD screen.* Gene hits (FDR \< 0.05) within one of the six highlighted gene groups are listed in the legend. Genes are arranged alphabetically along the horizontal axis. The lower right panel shows reporter mRNA level of the highlighted hits for the PTC reporter relative to the control mCherry-puro reporter. Box plot indicates median (center), interquartile range (box), and 3x interquartile range (whiskers) of the relative reporter mRNA level across all detected genes in the experiment. **c.** *Chemical modifier screen.* Cells were treated with 200 nM ISRIB or DMSO for 48 hours after 5 days of Cas9 induction. **d.** *Genetic modifier screen with GCN2 depletion.* Cells were transduced with lentiviruses expressing either a GCN2-targeting sgRNA or a control sgRNA, followed by Cas9 induction for 7 days. mRNA fold-changes and hits in c and d were calculated using MaGeCK as in b but with an FDR threshold of 0.01. Vertical axes in b, c and d indicate P-values from a permutation test as calculated by MAGeCK. Marker colors and shapes in c and d denote the highlighted gene groups from b.

### Figure 5

**GCN1 regulates cellular responses to an anti-leukemic drug.** **a.** *Chemogenomic ReLiC screen using homoharringtonine (HHT).* ReLiC-RBP cell pool with an EYFP reporter was treated with 1 *μ*M HHT or DMSO for 6 hours after 7 days of Cas9 induction. **b.** *GCN1 regulates mRNA levels upon HHT treatment.* Each point represents a gene knockout. Vertical axis indicates P-values from a permutation test as calculated by MAGeCK. **c.** *mRNA level changes upon HHT treatment for factors known to resolve ribosome collisions.* Each point represents a distinct sgRNA pair from the ReLiC screen shown in b. P-values comparing the indicated perturbations to cells expressing the nontargeting NLuc control sgRNA are from a two-sided t-test: \*\* (0.001 \< P \< 0.01), ns (P \> 0.05). **d.** *Immunoblots for phosphorylation of p38 in HEK293T cells +/- GCN1.* Cells were treated with HHT (1 *μ*M), anisomycin (ANS, 10 *μ*M), or DMSO for 1 hour. **e.** *GCN1-dependent changes in endogenous mRNA expression upon HHT treatment.* RNA-Seq was performed after treatment with HHT (1 *μ*M) or DMSO for 6 hours. Each points corresponds to a gene and represents the ratio of its mRNA levels between HHT and DMSO treatment. Black highlighted points correspond to immediate early genes (IEGs), which are also shown separately in the lower panel. **f.** *Metagene alignment of ribosome P-site density in 5<sup>′</sup> UTR and CDS across all detected transcripts.* Ribosome profiling was performed on +/- GCN1 cells after HHT (1*μ*M) or DMSO treatment for 1 hour. **g.** *Ribosome P-site density in 5<sup>′</sup> UTR and CDS of JUN and MYC transcripts.* Horizontal axis indicates position along the transcript in nucleotides.

### Extended Data Figure 1

**ReLiC library design and validation.**  
**a.**  *Validation of Cas9 activity in U2OS.*
sgEYFP and sgCTRL are single guide RNAs targeting EYFP or a non-targeting control, respectively.
Each histogram represents fluorescence of 10,000 cells as measured by flow cytometry.
`Days post Cas9' refers to days after addition of doxycycline to induce Cas9 expression.  
**b.**  *Comparison of integration into 293T and U2OS landing pads.*
BFP and mCherry fluorescence were measured for 10,000 cells, depicted as individual points.
Proportion of cells that are mCherry+ and BFP- (orange points) is indicated.
No cells in either parental control are mCherry+ and BFP-.  
**c.** *Depiction of cloning scheme for ReLiC library and reporters.*  
**d.** *Distribution of barcode read counts for sgRNA pairs in mRNA and genomic DNA.*  
**e.** *Number of unique barcodes linked to each sgRNA in ReLiC library.*  
**f.** *Correlation between distinct barcode sets in ReLiC fitness screens.*
Each point represents a unique sgRNA pair from the ReLiC RBP library.
For each sgRNA pair, individual linked barcodes were randomly partitioned into two sets of equal size (or to within a barcode for odd number of detected barcodes).
*r* refers to Pearson correlation coefficient between the barcode sets.

### Extended Data Figure 2

**Polysome ReLiC screen for regulators of mRNA translation.**  
**a.** *Correlation between replicates.*
Points represent individual sgRNAs in the ReLiC library.
Polysome to monosome ratios are median-centered across sgRNAs in the library.
*r* refers to Pearson correlation coefficient.  
**b.** *Change in polysome to supernatant ratio for ribosomal protein and ribosome biogenesis genes.*
Closed circles are genes that we call as gene hits (FDR < 0.05 with 3 or more concordant sgRNAs).
Open circles are genes that do not pass our gene hit threshold.  
**c.** *Gene ontology analysis of perturbations that decrease heavy polysome to monosome ratio.*  
Gene ontology analysis performed using GOrilla [@Eden2009] and a subset of enriched terms representative of specific gene classes are shown.  
**d.** *Comparison of heavy polysome to monosome and heavy polysome to supernatant ratios for selected translation-related factors.*  
**e.** *Polysome profiles of cell lines depleted of screen hits*
Profiles are normalized by 80S peak height.
P/M indicates ratio of area under the curve for polysome fractions to monosome fractions.  
**f.** *Comparison of heavy polysome to monosome ratio with growth fitness measured by mRNA and genomic DNA barcode seqencing 13 days after Cas9 induction for all gene knockouts.*  
**g.** *Comparison of heavy polysome to monosome ratio with growth fitness measured by genomic DNA barcode sequencing for gene knockouts in specific groups.*
Points correspond to genes targeted in the ReLiC-RBP library.
Shaded areas correspond to 95% confidence intervals for a linear fit of polysome to monosome ratio to growth fitness within each gene group.

### Extended Data Figure 3

**Isoform-specific splicing screen using ReLiC.**  
**a.** *Number of gene hits that increase the level of the indicated reporter isoform on indicated days after Cas9 induction.*  
**b.** *Correlation between barcode sets.*
For each sgRNA, individual linked barcodes were randomly partitioned into two sets, as in Extended Data Fig. 1d.
Each point represents a unique gene that was classified as a hit either with barcode Set A or barcode set B.
*r* refers to Pearson correlation coefficient between barcode sets.  
**c.** *Correlation between relative levels of different mRNA isoforms.*
Values represent Pearson correlation coefficients for pairwise comparison between the two barcode sets in B.  
**d.** *Depletion of genomic DNA barcodes corresponding to SF3b complex subunits after Cas9 induction.*

### Extended Data Figure 4

**Dissecting mRNA quality control using ReLiC.**  
**a.** *Validation of β-globin NMD reporters.* 
Relative reporter mRNA levels measured by qPCR (n=3).
Y-axis represents -ΔΔC~t~ value of indicated reporter mRNA relative to mCherry-Puro control mRNA.  
**b.** *Gene ontology analysis of perturbations that increase PTC reporter mRNA levels.*    
**c.** *Volcano plot of reporter mRNA levels with dual barcode screen.*  
Each point corresponds to a gene targeted by the ReLiC library.
Marker shape and color denotes one of highlighted gene groups.
Genes with FDR < 0.05 and belonging to one of the highlighted groups are listed in the legend.   
**d.** *PTC reporter levels for individual translation initiation complex subunits.*
Points denote mean and error bars denote standard deviation across sgRNAs for each gene.
P-values are as calculated by MAGeCK.   
**e.** *Growth fitness after depletion of translation initiation complex subunits.*

### Extended Data Figure 5

**GCN1 regulates cellular responses to the anti-leukemic drug homoharringtonine.**  
**a.** *Regulation of EYFP reporter levels by GCN2 after HHT treatment.*
Cell lines were treated with 1 μM GCN2i for 30m prior to 6h of 1 μM HHT treatment.
Y-axis represents the ratio of EYFP reporter barcode counts during indicated treatment compared to the DMSO-treated control in cells expressing indicated sgRNA.
**b.** *ZAK-dependent phosphorylation of p38 in HEK293T cells +/- GCN1 treated with HHT.* 
Cells were treated with nilotinib (1 μM) or DMSO for 30 minutes prior to addition of homoharringtonine (1 μM) treatment or DMSO for 1 hour.  
**c.** *GCN1-dependent changes to endogenous mRNA expression after HHT treatment in U937 cells.*
U937 cell lines were treated with indicated HHT concentrations or DMSO as a vehicle control for 6h.
Y-axis represents -ΔΔC~t~ value of either *EGR1* or *JUN* mRNA relative to *GAPDH* mRNA as measured by RT-qPCR (n = 3).  
**d.** *Regulation of endogenous mRNA expression by GCN2 and ZAK after HHT treatment.*
Cell lines were treated with 1 μM GCN2i or 1 μM of the ZAK inhibitors nilotinib and vemurafenib for 30m prior to 6h of 1 μM HHT treatment.
Y-axis represents -ΔΔC~t~ value of either *EGR1* or *JUN* mRNA relative to *GAPDH* mRNA as measured by RT-qPCR (n=3).  
**e.** *Polysome profiles of GCN1-depleted and control cell lines after HHT treatment.*  
Cells were treated with 1 μM HHT or DMSO for 1 hour prior to lysis. 
Polysome lysates were digested with 1 U micrococcal nuclease / μg of RNA prior to sucrose gradient sedimentation to isolate RNAse-resistant monosomes and disomes.  
**f.** *Ribosome P-site density on JUN and MYC mRNAs from previous ribosome profiling studies using harringtonine or lactimidomycin to arrest initiating ribosomes.*

## References

<div id="refs" class="references csl-bib-body" line-spacing="2">

<div id="ref-Gerstberger2014" class="csl-entry">

<span class="csl-left-margin">1. </span><span class="csl-right-inline">Gerstberger, S., Hafner, M. & Tuschl, T. [A census of human RNA-binding proteins](https://doi.org/10.1038/nrg3813). *Nature Reviews. Genetics* **15**, 829–845 (2014).</span>

</div>

<div id="ref-Hentze2018" class="csl-entry">

<span class="csl-left-margin">2. </span><span class="csl-right-inline">Hentze, M. W., Castello, A., Schwarzl, T. & Preiss, T. [A brave new world of RNA-binding proteins](https://doi.org/10.1038/nrm.2017.130). *Nature Reviews Molecular Cell Biology* **19**, 327–341 (2018).</span>

</div>

<div id="ref-Gebauer2021" class="csl-entry">

<span class="csl-left-margin">3. </span><span class="csl-right-inline">Gebauer, F., Schwarzl, T., Valcárcel, J. & Hentze, M. W. [RNA-binding proteins in human genetic disease](https://doi.org/10.1038/s41576-020-00302-y). *Nature Reviews Genetics* **22**, 185–198 (2021).</span>

</div>

<div id="ref-VanNostrand2020" class="csl-entry">

<span class="csl-left-margin">4. </span><span class="csl-right-inline">Van Nostrand, E. L. *et al.* [A large-scale binding and functional map of human RNA-binding proteins](https://doi.org/10.1038/s41586-020-2077-3). *Nature* **583**, 711–719 (2020).</span>

</div>

<div id="ref-Ray2013" class="csl-entry">

<span class="csl-left-margin">5. </span><span class="csl-right-inline">Ray, D. *et al.* [A compendium of RNA-binding motifs for decoding gene regulation](https://doi.org/10.1038/nature12311). *Nature* **499**, 172–177 (2013).</span>

</div>

<div id="ref-Schneider-lunitz2021" class="csl-entry">

<span class="csl-left-margin">6. </span><span class="csl-right-inline">Schneider-Lunitz, V., Ruiz-Orera, J., Hubner, N. & Heesch, S. van. [Multifunctional RNA-binding proteins influence <span class="nocase">mRNA</span> abundance and translational efficiency of distinct sets of target genes](https://doi.org/10.1371/journal.pcbi.1009658). *PLoS Computational Biology* **17**, e1009658 (2021).</span>

</div>

<div id="ref-England2022" class="csl-entry">

<span class="csl-left-margin">7. </span><span class="csl-right-inline">England, W. E. *et al.* [An atlas of posttranslational modifications on RNA binding proteins](https://doi.org/10.1093/nar/gkac243). *Nucleic Acids Research* **50**, 4329–4339 (2022).</span>

</div>

<div id="ref-Przybyla2022" class="csl-entry">

<span class="csl-left-margin">8. </span><span class="csl-right-inline">Przybyla, L. & Gilbert, L. A. [A new era in functional genomics screens](https://doi.org/10.1038/s41576-021-00409-w). *Nature Reviews Genetics* **23**, 89–103 (2022).</span>

</div>

<div id="ref-Genolet2022" class="csl-entry">

<span class="csl-left-margin">9. </span><span class="csl-right-inline">Genolet, O., Ravid Lustig, L. & Schulz, E. G. [Dissecting Molecular Phenotypes Through FACS-Based Pooled CRISPR Screens](https://doi.org/10.1007/7651_2021_457). *Methods in Molecular Biology (Clifton, N.J.)* **2520**, 1–24 (2022).</span>

</div>

<div id="ref-Gonatopoulos-Pournatzis2018" class="csl-entry">

<span class="csl-left-margin">10. </span><span class="csl-right-inline">Gonatopoulos-Pournatzis, T. *et al.* [Genome-wide CRISPR-Cas9 Interrogation of Splicing Networks Reveals a Mechanism for Recognition of Autism-Misregulated Neuronal Microexons](https://doi.org/10.1016/j.molcel.2018.10.008). *Molecular Cell* **72**, 510–524.e12 (2018).</span>

</div>

<div id="ref-Scarborough2021" class="csl-entry">

<span class="csl-left-margin">11. </span><span class="csl-right-inline">Scarborough, A. M. *et al.* [SAM homeostasis is regulated by CFIm-mediated splicing of MAT2A](https://doi.org/10.7554/eLife.64930). *eLife* **10**, e64930 (2021).</span>

</div>

<div id="ref-El-Brolosy2019" class="csl-entry">

<span class="csl-left-margin">12. </span><span class="csl-right-inline">El-Brolosy, M. A. *et al.* [Genetic compensation triggered by mutant <span class="nocase">mRNA</span> degradation](https://doi.org/10.1038/s41586-019-1064-z). *Nature* **568**, 193–197 (2019).</span>

</div>

<div id="ref-Replogle2022" class="csl-entry">

<span class="csl-left-margin">13. </span><span class="csl-right-inline">Replogle, J. M. *et al.* [Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq](https://doi.org/10.1016/j.cell.2022.05.013). *Cell* **185**, 2559–2575.e28 (2022).</span>

</div>

<div id="ref-Kowalski2024" class="csl-entry">

<span class="csl-left-margin">14. </span><span class="csl-right-inline">Kowalski, M. H. *et al.* [Multiplexed single-cell characterization of alternative polyadenylation regulators](https://doi.org/10.1016/j.cell.2024.06.005). *Cell* **187**, 4408–4425.e23 (2024).</span>

</div>

<div id="ref-Xu2024" class="csl-entry">

<span class="csl-left-margin">15. </span><span class="csl-right-inline">Xu, Z., Sziraki, A., Lee, J., Zhou, W. & Cao, J. [Dissecting key regulators of transcriptome kinetics through scalable single-cell RNA profiling of pooled CRISPR screens](https://doi.org/10.1038/s41587-023-01948-9). *Nature Biotechnology* **42**, 1218–1223 (2024).</span>

</div>

<div id="ref-Herholt2018" class="csl-entry">

<span class="csl-left-margin">16. </span><span class="csl-right-inline">Herholt, A. *et al.* [Pathway sensor-based functional genomics screening identifies modulators of neuronal activity](https://doi.org/10.1038/s41598-018-36008-9). *Scientific Reports* **8**, 17597 (2018).</span>

</div>

<div id="ref-Muller2020" class="csl-entry">

<span class="csl-left-margin">17. </span><span class="csl-right-inline">Muller, R., Meacham, Z. A., Ferguson, L. & Ingolia, N. T. [CiBER-seq dissects genetic networks by quantitative CRISPRi profiling of expression phenotypes](https://doi.org/10.1126/science.abb9662). *Science* **370**, (2020).</span>

</div>

<div id="ref-Alford2021" class="csl-entry">

<span class="csl-left-margin">18. </span><span class="csl-right-inline">Alford, B. D. *et al.* [ReporterSeq reveals genome-wide dynamic modulators of the heat shock response across diverse stressors](https://doi.org/10.7554/eLife.57376). *eLife* **10**, e57376 (2021).</span>

</div>

<div id="ref-Sack2016" class="csl-entry">

<span class="csl-left-margin">19. </span><span class="csl-right-inline">Sack, L. M., Davoli, T., Xu, Q., Li, M. Z. & Elledge, S. J. [Sources of Error in Mammalian Genetic Screens](https://doi.org/10.1534/g3.116.030973). *G3 (Bethesda, Md.)* **6**, 2781–2790 (2016).</span>

</div>

<div id="ref-Ellis2005" class="csl-entry">

<span class="csl-left-margin">20. </span><span class="csl-right-inline">Ellis, J. [Silencing and variegation of gammaretrovirus and lentivirus vectors](https://doi.org/10.1089/hum.2005.16.1241). *Human Gene Therapy* **16**, 1241–1246 (2005).</span>

</div>

<div id="ref-Warner1963" class="csl-entry">

<span class="csl-left-margin">21. </span><span class="csl-right-inline">Warner, J. R., Knopf, P. M. & Rich, A. [A multiple ribosomal structure in protein synthesis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC300639). *Proceedings of the National Academy of Sciences of the United States of America* **49**, 122–129 (1963).</span>

</div>

<div id="ref-Noll1963" class="csl-entry">

<span class="csl-left-margin">22. </span><span class="csl-right-inline">Noll, H., Staehelin, T. & Wettstein, F. O. [Ribosomal Aggregates Engaged in Protein Synthesis : Ergosome Breakdown and Messenger Ribonucleic Acid Transport](https://doi.org/10.1038/198632a0). *Nature* **198**, 632–638 (1963).</span>

</div>

<div id="ref-Gierer1963" class="csl-entry">

<span class="csl-left-margin">23. </span><span class="csl-right-inline">Gierer, A. [Function of aggregated reticulocyte ribosomes in protein synthesis](https://doi.org/10.1016/S0022-2836(63)80131-X). *Journal of Molecular Biology* **6**, 148–IN8 (1963).</span>

</div>

<div id="ref-Zhang1998" class="csl-entry">

<span class="csl-left-margin">24. </span><span class="csl-right-inline">Zhang, J., Sun, X., Qian, Y. & Maquat, L. E. [Intron function in the nonsense-mediated decay of beta-globin <span class="nocase">mRNA</span>: Indications that pre-<span class="nocase">mRNA</span> splicing in the nucleus can influence <span class="nocase">mRNA</span> translation in the cytoplasm.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1369660) *RNA* **4**, 801–815 (1998).</span>

</div>

<div id="ref-Li2014" class="csl-entry">

<span class="csl-left-margin">25. </span><span class="csl-right-inline">Li, W. *et al.* [MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens](https://doi.org/10.1186/s13059-014-0554-4). *Genome Biology* **15**, 554 (2014).</span>

</div>

<div id="ref-Wagner2016" class="csl-entry">

<span class="csl-left-margin">26. </span><span class="csl-right-inline">Wagner, S., Herrmannová, A., Šikrová, D. & Valášek, L. S. [Human <span class="nocase">eIF3b</span> and <span class="nocase">eIF3a</span> serve as the nucleation core for the assembly of <span class="nocase">eIF3</span> into two interconnected modules: The yeast-like core and the octamer](https://doi.org/10.1093/nar/gkw972). *Nucleic Acids Research* **44**, 10772–10788 (2016).</span>

</div>

<div id="ref-Darnell2018" class="csl-entry">

<span class="csl-left-margin">27. </span><span class="csl-right-inline">Darnell, A. M., Subramaniam, A. R. & O’Shea, E. K. [Translational Control through Differential Ribosome Pausing during Amino Acid Limitation in Mammalian Cells](https://doi.org/10.1016/j.molcel.2018.06.041). *Molecular Cell* **71**, 229–243.e11 (2018).</span>

</div>

<div id="ref-Collart2016" class="csl-entry">

<span class="csl-left-margin">28. </span><span class="csl-right-inline">Collart, M. A. [The Ccr4‐Not complex is a key regulator of eukaryotic gene expression](https://doi.org/10.1002/wrna.1332). *WIREs RNA* **7**, 438–454 (2016).</span>

</div>

<div id="ref-Harper2016" class="csl-entry">

<span class="csl-left-margin">29. </span><span class="csl-right-inline">Harper, J. W. & Bennett, E. J. [Proteome complexity and the forces that drive proteome imbalance](https://doi.org/10.1038/nature19947). *Nature* **537**, 328–338 (2016).</span>

</div>

<div id="ref-Frydman2012" class="csl-entry">

<span class="csl-left-margin">30. </span><span class="csl-right-inline">Frydman, J. [Mechanism and Function of the Eukaryotic Ring-Shaped Chaperonin TRiC/CCT](https://doi.org/10.1016/j.bpj.2011.11.2340). *Biophysical Journal* **102**, 428a (2012).</span>

</div>

<div id="ref-Balakrishnan2022" class="csl-entry">

<span class="csl-left-margin">31. </span><span class="csl-right-inline">Balakrishnan, R. *et al.* [Principles of gene regulation quantitatively connect DNA to RNA and proteins in bacteria](https://doi.org/10.1126/science.abk2066). *Science (New York, N.Y.)* **378**, eabk2066 (2022).</span>

</div>

<div id="ref-Metzl-Raz2017" class="csl-entry">

<span class="csl-left-margin">32. </span><span class="csl-right-inline">Metzl-Raz, E. *et al.* [Principles of cellular resource allocation revealed by condition-dependent proteome profiling](https://doi.org/10.7554/eLife.28034). *eLife* **6**, e28034 (2017).</span>

</div>

<div id="ref-Juszkiewicz2020" class="csl-entry">

<span class="csl-left-margin">33. </span><span class="csl-right-inline">Juszkiewicz, S., Speldewinde, S. H., Wan, L., Svejstrup, J. Q. & Hegde, R. S. [The ASC-1 Complex Disassembles Collided Ribosomes](https://doi.org/10.1016/j.molcel.2020.06.006). *Molecular Cell* **79**, 603–614.e8 (2020).</span>

</div>

<div id="ref-Datta1988" class="csl-entry">

<span class="csl-left-margin">34. </span><span class="csl-right-inline">Datta, B., Chakrabarti, D., Roy, A. L. & Gupta, N. K. [Roles of a 67-<span class="nocase">kDa</span> polypeptide in reversal of protein synthesis inhibition in heme-deficient reticulocyte lysate](https://doi.org/10.1073/pnas.85.10.3324). *Proceedings of the National Academy of Sciences of the United States of America* **85**, 3324–3328 (1988).</span>

</div>

<div id="ref-Wahl2015" class="csl-entry">

<span class="csl-left-margin">35. </span><span class="csl-right-inline">Wahl, M. C. & Lührmann, R. [SnapShot: Spliceosome Dynamics I](https://doi.org/10.1016/j.cell.2015.05.050). *Cell* **161**, 1474–e1 (2015).</span>

</div>

<div id="ref-Hluchy2022" class="csl-entry">

<span class="csl-left-margin">36. </span><span class="csl-right-inline">Hluchý, M. *et al.* [CDK11 regulates pre-<span class="nocase">mRNA</span> splicing by phosphorylation of SF3B1](https://doi.org/10.1038/s41586-022-05204-z). *Nature* **609**, 829–834 (2022).</span>

</div>

<div id="ref-Cheng1994" class="csl-entry">

<span class="csl-left-margin">37. </span><span class="csl-right-inline">Cheng, J., Belgrader, P., Zhou, X. & Maquat, L. E. [Introns are cis Effectors of the Nonsense-Codon-Mediated Reduction in Nuclear <span class="nocase">mRNA</span> Abundance](https://doi.org/10.1128/mcb.14.9.6317-6325.1994). *Molecular and Cellular Biology* **14**, 6317–6325 (1994).</span>

</div>

<div id="ref-Darman2015" class="csl-entry">

<span class="csl-left-margin">38. </span><span class="csl-right-inline">Darman, R. B. *et al.* [Cancer-Associated SF3B1 Hotspot Mutations Induce Cryptic 3’ Splice Site Selection through Use of a Different Branch Point](https://doi.org/10.1016/j.celrep.2015.09.053). *Cell Reports* **13**, 1033–1045 (2015).</span>

</div>

<div id="ref-Cretu2016" class="csl-entry">

<span class="csl-left-margin">39. </span><span class="csl-right-inline">Cretu, C. *et al.* [Molecular Architecture of SF3b and Structural Consequences of Its Cancer-Related Mutations](https://doi.org/10.1016/j.molcel.2016.08.036). *Molecular Cell* **64**, 307–319 (2016).</span>

</div>

<div id="ref-Wang2011" class="csl-entry">

<span class="csl-left-margin">40. </span><span class="csl-right-inline">Wang, D. *et al.* [Inhibition of Nonsense-Mediated RNA Decay by the Tumor Microenvironment Promotes Tumorigenesis](https://doi.org/10.1128/MCB.05704-11). *Molecular and Cellular Biology* **31**, 3670–3680 (2011).</span>

</div>

<div id="ref-Sidrauski2013" class="csl-entry">

<span class="csl-left-margin">41. </span><span class="csl-right-inline">Sidrauski, C. *et al.* [Pharmacological brake-release of <span class="nocase">mRNA</span> translation enhances cognitive memory](https://doi.org/10.7554/eLife.00498). *eLife* **2**, e00498 (2013).</span>

</div>

<div id="ref-Wek1995" class="csl-entry">

<span class="csl-left-margin">42. </span><span class="csl-right-inline">Wek, S. A., Zhu, S. & Wek, R. C. [The histidyl-<span class="nocase">tRNA</span> synthetase-related sequence in the <span class="nocase">eIF</span>-2 alpha protein kinase GCN2 interacts with <span class="nocase">tRNA</span> and is required for activation in response to starvation for different amino acids](https://doi.org/10.1128/MCB.15.8.4497). *Molecular and Cellular Biology* **15**, 4497–4506 (1995).</span>

</div>

<div id="ref-Dong2000" class="csl-entry">

<span class="csl-left-margin">43. </span><span class="csl-right-inline">Dong, J., Qiu, H., Garcia-Barrio, M., Anderson, J. & Hinnebusch, A. G. [Uncharged <span class="nocase">tRNA</span> activates GCN2 by displacing the protein kinase moiety from a bipartite <span class="nocase">tRNA</span>-binding domain](https://doi.org/10.1016/s1097-2765(00)00028-9). *Molecular Cell* **6**, 269–279 (2000).</span>

</div>

<div id="ref-Gandhi2014" class="csl-entry">

<span class="csl-left-margin">44. </span><span class="csl-right-inline">Gandhi, V., Plunkett, W. & Cortes, J. E. [Omacetaxine: A protein translation inhibitor for treatment of chronic myelogenous leukemia](https://doi.org/10.1158/1078-0432.CCR-13-1283). *Clinical Cancer Research: An Official Journal of the American Association for Cancer Research* **20**, 1735–1740 (2014).</span>

</div>

<div id="ref-Fresno1977" class="csl-entry">

<span class="csl-left-margin">45. </span><span class="csl-right-inline">Fresno, M., Jiménez, A. & Vázquez, D. [Inhibition of translation in eukaryotic systems by harringtonine](https://doi.org/10.1111/j.1432-1033.1977.tb11256.x). *European Journal of Biochemistry* **72**, 323–330 (1977).</span>

</div>

<div id="ref-Ingolia2011" class="csl-entry">

<span class="csl-left-margin">46. </span><span class="csl-right-inline">Ingolia, N. T., Lareau, L. F. & Weissman, J. S. [Ribosome Profiling of Mouse Embryonic Stem Cells Reveals the Complexity and Dynamics of Mammalian Proteomes](https://doi.org/10.1016/j.cell.2011.10.002). *Cell* **147**, 789–802 (2011).</span>

</div>

<div id="ref-Marton1993" class="csl-entry">

<span class="csl-left-margin">47. </span><span class="csl-right-inline">Marton, M. J., Crouch, D. & Hinnebusch, A. G. [GCN1, a translational activator of GCN4 in Saccharomyces cerevisiae, is required for phosphorylation of eukaryotic translation initiation factor 2 by protein kinase GCN2](https://doi.org/10.1128/mcb.13.6.3541-3556.1993). *Molecular and Cellular Biology* **13**, 3541–3556 (1993).</span>

</div>

<div id="ref-Pochopien2021" class="csl-entry">

<span class="csl-left-margin">48. </span><span class="csl-right-inline">Pochopien, A. A. *et al.* [Structure of Gcn1 bound to stalled and colliding 80S ribosomes](https://doi.org/10.1073/pnas.2022756118). *Proceedings of the National Academy of Sciences* **118**, e2022756118 (2021).</span>

</div>

<div id="ref-Oltion2023" class="csl-entry">

<span class="csl-left-margin">49. </span><span class="csl-right-inline">Oltion, K. *et al.* [An E3 ligase network engages GCN1 to promote degradation of translation factors on stalled ribosomes](https://doi.org/10.1016/j.cell.2022.12.025). *Cell* **186**, 346–362.e17 (2023).</span>

</div>

<div id="ref-Muller2023" class="csl-entry">

<span class="csl-left-margin">50. </span><span class="csl-right-inline">Müller, M. B. D., Kasturi, P., Jayaraj, G. G. & Hartl, F. U. [Mechanisms of readthrough mitigation reveal principles of GCN1-mediated translational quality control](https://doi.org/10.1016/j.cell.2023.05.035). *Cell* **186**, 3227–3244.e20 (2023).</span>

</div>

<div id="ref-Wu2020" class="csl-entry">

<span class="csl-left-margin">51. </span><span class="csl-right-inline">Wu, C. C.-C., Peterson, A., Zinshteyn, B., Regot, S. & Green, R. [Ribosome collisions trigger general stress responses to regulate cell fate](https://doi.org/10.1016/j.cell.2020.06.006). *Cell* **182**, 404–416.e14 (2020).</span>

</div>

<div id="ref-Sinha2020" class="csl-entry">

<span class="csl-left-margin">52. </span><span class="csl-right-inline">Sinha, N. K. *et al.* [EDF1 coordinates cellular responses to ribosome collisions](https://doi.org/10.7554/eLife.58828). *eLife* **9**, e58828 (2020).</span>

</div>

<div id="ref-Hann1988" class="csl-entry">

<span class="csl-left-margin">53. </span><span class="csl-right-inline">Hann, S. R., King, M. W., Bentley, D. L., Anderson, C. W. & Eisenman, R. N. [A non-AUG translational initiation in c-myc exon 1 generates an N-terminally distinct protein whose synthesis is disrupted in Burkitt’s lymphomas](https://doi.org/10.1016/0092-8674(88)90507-7). *Cell* **52**, 185–195 (1988).</span>

</div>

<div id="ref-Short2002" class="csl-entry">

<span class="csl-left-margin">54. </span><span class="csl-right-inline">Short, J. D. & Pfarr, C. M. [Translational Regulation of the JunD Messenger RNA \*](https://doi.org/10.1074/jbc.M204553200). *Journal of Biological Chemistry* **277**, 32697–32705 (2002).</span>

</div>

<div id="ref-Gonzalez-Sanchez2024" class="csl-entry">

<span class="csl-left-margin">55. </span><span class="csl-right-inline">González-Sánchez, A. M., Castellanos-Silva, E. A., Díaz-Figueroa, G. & Cate, J. H. D. [JUN <span class="nocase">mRNA</span> translation regulation is mediated by multiple 5’ UTR and start codon features](https://doi.org/10.1371/journal.pone.0299779). *PLOS ONE* **19**, e0299779 (2024).</span>

</div>

<div id="ref-Natsume2016" class="csl-entry">

<span class="csl-left-margin">56. </span><span class="csl-right-inline">Natsume, T., Kiyomitsu, T., Saga, Y. & Kanemaki, M. T. [Rapid Protein Depletion in Human Cells by Auxin-Inducible Degron Tagging with Short Homology Donors](https://doi.org/10.1016/j.celrep.2016.03.001). *Cell Reports* **15**, 210–218 (2016).</span>

</div>

<div id="ref-Ingolia2012" class="csl-entry">

<span class="csl-left-margin">57. </span><span class="csl-right-inline">Ingolia, N. T., Brar, G. A., Rouskin, S., McGeachy, A. M. & Weissman, J. S. [The ribosome profiling strategy for monitoring translation in vivo by deep sequencing of ribosome-protected <span class="nocase">mRNA</span> fragments](https://doi.org/10.1038/nprot.2012.086). *Nature Protocols* **7**, 1534–1550 (2012).</span>

</div>

<div id="ref-Ingolia2009" class="csl-entry">

<span class="csl-left-margin">58. </span><span class="csl-right-inline">Ingolia, N. T., Ghaemmaghami, S., Newman, J. R. S. & Weissman, J. S. [Genome-Wide Analysis in Vivo of Translation with Nucleotide Resolution Using Ribosome Profiling](https://doi.org/10.1126/science.1168978). *Science* **324**, 218–223 (2009).</span>

</div>

<div id="ref-Lee2012" class="csl-entry">

<span class="csl-left-margin">59. </span><span class="csl-right-inline">Lee, S. *et al.* [Global mapping of translation initiation sites in mammalian cells at single-nucleotide resolution](https://doi.org/10.1073/pnas.1207846109). *Proceedings of the National Academy of Sciences* **109**, E2424–E2432 (2012).</span>

</div>

<div id="ref-Koster2012" class="csl-entry">

<span class="csl-left-margin">60. </span><span class="csl-right-inline">Köster, J. & Rahmann, S. [Snakemake—a scalable bioinformatics workflow engine](https://doi.org/10.1093/bioinformatics/bts480). *Bioinformatics* **28**, 2520–2522 (2012).</span>

</div>

</div>

# Methods

## Plasmid construction

Plasmids used in this study are listed in Supplementary Table 2, which also includes accession numbers for plasmids deposited to Addgene. Oligonucleotides used in this study are listed in Supplementary Table 3. Detailed cloning steps for all plasmid vectors constructed for this study are described in Supplementary Methods. DNA fragments used for cloning were either excised out by restriction digestion or amplified by PCR from suitable templates. Fragments were assembled together using isothermal assembly, and transformed into NEB10beta cells. All constructs were verified by restriction digestion and Sanger or long read sequencing.

## Cell culture

Cell lines used in this study are listed in Supplementary Table 4. HEK293T (RRID:CVCL_0063, ATCC CRL-3216) and U2OS (RRID:CVCL_0042, ATCC HTB-96) cells were cultured in Dulbecco’s modified Eagle medium (DMEM 1X, with 4.5 g/L D-glucose, + L-glutamine, - sodium pyruvate, Gibco 11965-092) supplemented with 10% FBS (Thermo 26140079) and 1X penicillin–streptomycin (Gibco 15140-122), and were passaged using 0.25% trypsin in EDTA (Gibco 25200-056). HEK293T cell lines were authenticated by short tandem repeat analysis. U2OS cells were not authenticated. U937 cells (RRID:CVCL_0007) were cultured in RPMI supplemented with 10% heat-inactivated FBS, penicillin-streptomycin, and 1X Glutamax (Gibco 35050-061). Cells were grown at 37C in 5% CO2. Cell lines were periodically confirmed to be free of mycoplasma contamination.

## Generation of landing pad cell lines

To generate an initial *attP* landing pad line, HEK293T cells were transfected with landing pad plasmid (pHPHS232) and pASHS29 (AAVS1 T2 CRISPR in pX330<sup>[56](#ref-Natsume2016)</sup>/Addgene 72833) using polyethylenimine. U2OS cells were nucleofected with the same plasmids using the Nucleofector V kit (Lonza) in a Nucleofector 2b system. HEK293T cells were selected with 10 *μ*g/ml Blasticidin S, added 96 hours post-transfection and maintained for 4 days. Blasticidin selection was omitted for U2OS cells. At this point, BFP expression was induced in all lines by adding 2 *μ*g/ml doxycycline. 24 hours after doxycycline induction, the cultures were enriched for BFP+ cells using a FACSAria II flow cytometer (BD). The sorted BFP+ U2OS cells were kept as a polyclonal cell line (hsPN279). HEK293T landing pad clones were isolated by limiting dilution into 96-well plates. After isolating clones, two were pooled into a single cell line (hsPB126).

To integrate a *Cas9* expression cassette with an orthogonal *attP\** site into the initial *attP* landing pad clonal lines, hsPB126 was transfected with Cas9 landing pad plasmid (pHPHS800) and Bxb1 expression plasmid (pHPHS115) using TransIT-LT1 reagent (Mirus) while hsPN279 was nucleofected with the same plasmids using the Nucleofector V kit. 72 hours post-transfection, *hygromycin phosphotransferase* was induced by adding 2 *μ*g/ml doxycycline, then cells were selected with 150 *μ*g/ml Hygromycin B, added 96 hours post-transfection. After 7 days, doxycycline and Hygromycin B were removed from cells. At this point, the HEK293T cells were further selected using 10 *μ*g/ml Blasticidin for 7 days, and this polyclonal cell line (hsPN266) was used for subsequent experiments. Instead of selecting with Blasticidin, the U2OS cells were further enriched for BFP+ cells using a FACSAria II flow cytometer (BD). The sorted BFP+ U2OS cells were kept as a polyclonal cell line (hsPN280).

## Validation of reporter integration and Cas9 activity in landing pad cell lines

hsPN266 (HEK293T *attP\* Cas9*) cells were seeded to 70% confluency in a 6-well or 12-well culture dish, and transfected with 1 *μ*g (6-well) or 200 ng (12-well) of Bxb1 expression vector (pHPHS115) and 2 *μ*g (6-well) or 800 ng (12-well) of *attB\**-containing sgFLUC-EYFP or sgYFP-EYFP reporter plasmid (pHPHS885 and pHPHS886, respectively) using Fugene HD (6-well) or Lipofectamine 2000 (12-well) reagent. 1x10<sup>6</sup> hsPN280 (U2OS *attP\* Cas9*) cells were nucleofected with 1 *μ*g of pHPHS115 and 2 *μ*g of pHPHS885 or pHPHS886 using the Nucleofector V kit. Each culture was selected with 2 *μ*g/ml puromycin, added 72 hours post-transfection. Flow cytometry was performed on a small split of cells from each culture not selected with puromycin to measure integration efficiency 7 days after transfection (Extended Data Fig. 1b) . Puromycin selection was ended after 4 days on these cell lines (referred to as hsPN349,351,353,354). 24h after ending puromycin selection, 2 *μ*g/ml doxycycline was added to induce Cas9 expression. After inducing Cas9, reporter expression was monitored over time by flow cytometry on the indicated days (Fig. 1b, Extended Data Fig. 1a).

## sgRNA insert-barcode linkage sequencing

sgRNA insert-barcode linkages were determined at the step right after barcodes were added to the cloned sgRNA plasmid pool, prior to adding AmpR between the sgRNAs. A 422 bp amplicon containing both sgRNAs and 20xN barcodes was generated from 1.5 ng of pHPHS932 plasmid by 10 cycles of PCR using oKC196/oPN726 primers and Phusion polymerase (Thermo). This product cut out from a 1.5% agarose gel and cleaned using the Zymoclean Gel DNA Recovery Kit (Zymo). This sample was sequenced on an Illumina NextSeq 2000 using custom sequencing primers: oAS1701 for Read 1 (26 cycles), oKC186 for Index 1 (6 cycles), oAS1702 for Index 2 (20 cycles), and oKC185 for Read 2 (75 cycles).

## Integration of plasmid libraries into landing pad

hsPN266 (HEK293T *attP\* Cas9*) cells were seeded to 60% confluency in one 15 cm dish per library. 20 *μ*g of *attB\**-containing reporter library plasmid (pAS243, pAS244, pHPHS951) and 5 *μ*g of Bxb1 expression vector (pHPHS115) were transfected per 15 cm dish using TransIT-LT1 reagent (Mirus). Each library was transfected into a single 15 cm dish then expanded into four 15 cm dishes 48 hours post-transfection. Cells were selected with 2 *μ*g/ml puromycin, added 72 hours post-transfection. Puromycin selection was ended after 4 days, and library cell lines (referred to as hsPN305, hsPN306, hsPN285) were contracted back into a single 15 cm dish. 24h after ending puromycin selection, 2 *μ*g/ml doxycycline was added to induce Cas9 expression, and libraries were expanded into three 15 cm dishes – one each for RNA and gDNA harvests the next day plus a third for continued propagation. This splitting procedure was repeated every other day from the propagation dish, so harvests could be taking throughout the duration of the screen. At no point were cultures bottlenecked to fewer than 5x10<sup>6</sup> cells.

## Library genomic DNA extraction

For each harvest, reporter library genomic DNA was harvested from one 50% confluent 15 cm dish of cells stably expressing the ReLiC library. Genomic DNA was harvested using Quick-DNA Miniprep kit (Zymo), following the manufacturer’s instructions, with 2.5 ml of genomic DNA lysis buffer per 15 cm plate. 30 *μ*g of purified genomic DNA from each library sample was sheared into \~350 nucleotide length fragments by sonication for 10minutes on ice using a Diagenode Bioruptor. Sheared gDNA was then *in vitro* transcribed into RNA (denoted gRNA below and in analysis code) starting from the T7 promoter region in the insert cassette using the HiScribe T7 High Yield RNA Synthesis Kit (NEB). Transcribed gRNA was cleaned using the RNA Clean and Concentrator kit (Zymo).

## Library mRNA extraction

For each harvest, reporter library mRNA was harvested from one 50-75% confluent 15 cm dish of cells stably expressing the ReLiC library. Total RNA was harvested by using 3.5 ml of Trizol reagent (Thermo) to lyse cells directly on the plate, and then RNA was extracted from these lysates using the Direct-zol RNA Miniprep kit (Zymo) following the manufacturer’s protocol. polyA+ mRNA was extracted from total RNA using oligo dT25 magnetic beads (NEB). 30-50 *μ*g of total RNA was used as polyA selection input for total barcode counting libraries from each sample while 10-12 *μ*g was used as input for splicing or polysome fraction barcode counting libraries. 4 *μ*l of oligo dT25 beads were used per 1 *μ*g of total RNA input.

## mRNA and genomic DNA barcode sequencing

100-500 ng of polyA-selected mRNA or *in vitro* transcribed gRNA from each library was reverse transcribed into cDNA using SuperScript IV reverse transcriptase (Thermo) following the manufacturer’s protocol. For RT, we used a primer that binds downstream of the 20xN reporter barcode: either oPN777 for mRNA barcode 1, oPN731 for gRNA barcode 1, or oPN779 for mRNA barcode 2. oPN777 and oPN779 contain a 7 nt UMI. Libraries for sequencing total levels of barcode 1 or barcode 2 in each sample were performed in a single step. For both barcodes, a 100-200 *μ*l PCR was performed using Phusion polymerase (Thermo) for 20-25 cycles with cDNA template comprising 1/5th of the final volume, and oPN776 was used as a constant reverse primer that binds the Illumina P5 sequence present on oPN777 and oPN779. Indexed forward primers that bind a constant region upstream of each barcode were used to enable pooled sequencing of different samples (one of oPN730, oPN738, oPN809, oPN815-822, or oJY1-14 for Barcode 1 or one of oPN734, oPN739, or oPN823-825 for Barcode 2). All of these reactions generated a 181 bp amplicon that was cut out from a 2% agarose gel and cleaned using the Zymoclean Gel DNA Recovery Kit (Zymo).

For splicing screens, two rounds of PCR were performed. Round 1 was performed as a 50 *μ*l PCR for 30 cycles, again with cDNA template comprising 1/5th of the final volume and oPN776 as a constant reverse primer. The forward primer for Round 1 was chosen based on the measured splicing event: oPN841 for intron 1 retention, oPN789 for intron 2 retention, or oAS2029 for exon 2 skipping. These generate 532, 302, and 286 bp amplicons, respectively, which were cut out from a 2% agarose gel and cleaned using the Zymoclean Gel DNA Recovery Kit (Zymo), eluting in 15 *μ*l. Round 2 PCR was then essentially the same as the single-step PCR for total Barcode 1 sequencing, except reactions were 20 *μ*l, used 4 *μ*l of cleaned Round 1 product as template, and proceeded for 5 cycles.

Libraries were sequenced on an Illumina NextSeq 2000 using custom sequencing primers. Custom primers for Barcode 1 were oAS1701 for Read 1, oPN732 for Index 1, oPN775 for Index 2, and oPN731 for Read 2. Custom primers for Barcode 2 were oPN735 for Read 1, oPN737 for Index 1, oPN778 for Index 2, and oPN736 for Read 2. Read lengths varied between sequencing runs with 10% phiX spiked in.

## Polysome fractionation for ReLiC

After Cas9 induction, 293T cells expressing ReLiC libraries were passaged for 6 days. On day 6, lysates were prepared from each library at 30 percent confluency in a 15 cm dish. Cultures were treated with 100 *μ*g/ml cycloheximide for 5 minutes prior to harvest, then cells were trypsinized (including 100 *μ*g/ml cycloheximide) and pelleted at 300g for 5 min. Cell pellets were lysed on ice in 300 *μ*l of polysome lysis buffer (10 mM Tris-HCl pH 7.4 (Ambion), 132 mM NaCl (Ambion), 1.4 mM MgCl2 (Ambion), 19 mM DTT (Sigma), 142 *μ*g/ml cycloheximide (Sigma), 0.1% Triton X-100 (Fisher), 0.2% NP-40 (Pierce), 607 U/ml SUPERase-In RNase Inhibitor (Invitrogen)) with periodic vortex mixing. Lysates were clarified by centrifugation at 9,300g for 5 min and supernatants were transferred to fresh tubes. This total lysate was split into two parts: 50 *μ*l for total mRNA isolation, and 250 *μ*l for polysome profiling. For each sample, the 250 *μ*L lysate fraction was layered onto a 10%–50% (w/v) linear sucrose gradient (Fisher) containing 2 mM DTT (Sigma) and 100 *μ*g/mL heparin (Sigma). The gradients were centrifuged at 235,000g (37,000 rpm) for 2.5 h at 4°C in a Beckman SW41Ti rotor in Seton 7030 ultracentrifuge tubes. After centrifugation, samples were fractionated using a Biocomp Gradient Station by upward displacement into collection tubes, through a Bio-Rad EM-1 UV monitor (Bio-Rad) for continuous measurement of the absorbance at 260 nm. 820 *μ*l of TRIzol Reagent (Invitrogen) were added to each RNA fraction. Total (input), supernatant (fraction 1 and 2), monosome-associated (fraction 4 and 5), low polysome-associated (fractions 6-9), and high polysome-associated (fractions 10-13) samples were pooled as necessary and RNA isolated from TRIzol (Invitrogen) using the Direct-zol RNA Miniprep Plus Kit (Zymo Research) with DNaseI treatment according to manufacturer’s directions. These RNA samples were then subject to barcode sequencing as described above.

## Polysome profiling

To examine the effect of CNOT1, PSMA4, RPS19, TCP1, and FLUC depletion on global translation, polysome profiling was performed with the five cell lines from “CRISPR-Cas9 mediated gene knockout for polysome profiling” (Supplementary Methods). Polysome profiling was performed similar to the polysome ReLiC scree, except each sample was harvested from one 10-cm dish at \~60% confluency. To examine whether GCN1 affects the level of RNAse-resistant disomes during HHT treatment, polysome profiling was performed with four different samples: 293T cells expressing sgGCN1 and sgFLUC from “CRISPR-Cas9 mediated gene knockout for RNA-seq” (Supplementary Methods) after 1 week of Cas9 induction, treated for 1 hour with 1 *μ*M HHT or DMSO. Polysome profiling was performed similar to the individual gene peturbations from above, except prior to loading onto sucrose gradients, lysates were incubated with or without the addition of 1 U of micrococcal nuclease per *μ*g of RNA and 5 *μ*M CaCl<sub>2</sub> at room temperature for 1 hour. Micrococcal nuclease digests were quenched by addition of 5 *μ*M EGTA prior to loading on sucrose gradients.

## RNA-seq

Total RNA was isolated using the Direct-zol RNA Miniprep kit (Zymo). polyA+ mRNA was extracted from total RNA using oligo dT25 magnetic beads (NEB). 2 *μ*g of total RNA input and 25 *μ*l of oligo dT25 beads were used per sample. Sequencing libraries were generated from polyA+ mRNA using the NEBNext Ultra II Directional RNA Library Prep Kit (NEB) and sequenced on a NextSeq 2000 (Illumina) with 2x50 cycle paired-end reads.

## Ribosome profiling

Ribosome profiling was performed with four different samples: 293T cells expressing sgGCN1 and sgFLUC from “CRISPR-Cas9 mediated gene knockout for RNA-seq” (Supplementary Methods) after 1 week of Cas9 induction, treated for 1 hour with 1 *μ*M HHT or DMSO. For each sample, we used one 15-cm plate of cells, seeded to \~40% confluence at harvest. Ribosome profiling protocol was adapted from<sup>[57](#ref-Ingolia2012)</sup> with the following modifications. For sample harvesting, we removed media from each plate and flash froze samples by placing the plate in liquid nitrogen and transferred to -80C until lysis. We performed nuclease footprinting treatment by adding 80 U RNase I (Invitrogen AM2294) to 25 *μ*g of RNA. We gel-purified ribosome protected fragments with length between 26 and 34 nucleotides using RNA oligo size markers. We used polyA tailing instead of linker ligation following previous studies<sup>[58](#ref-Ingolia2009),[59](#ref-Lee2012)</sup>. Libraries were sequenced on an Illumina Nextseq 2000 in 50bp single end mode.

## Immunoblotting

sgGCN1- and sgFLUC-expressing cell lines hsPN309 and hsPN313 were passaged for 1 week with 2 *μ*g/ml doxycycline to deplete GCN1. These lines were subsequently incubated with homoharringtonine, anisomycin, or DMSO at indicated concentrations for 60 min before harvest. If indicated, cells were pre-treated with nilotinib for 30 mins prior to starting HHT treatment. Homoharringtonine, anisomycin, and nilotinib were dissolved in DMSO. Cells were rinsed with PBS and lysed in RIPA buffer. Lysates were kept on ice during preparation and clarified by centrifugation at 21,000g (15,000 rpm) for 10 min. After clarification, supernatants were boiled in Laemmli loading buffer containing DTT, and Western blots were performed using standard molecular biology procedures (see Supplementary Figure 1 and Supplementary Figure 2 for unprocessed blot images). Proteins were resolved by 4%–20% Criterion TGX protein gels (Bio-Rad) and transferred to PVDF membranes using a Trans-Blot Turbo transfer system (Bio-Rad). Membranes were blocked with 5% BSA (Thermo) in TBST and incubated with primary antibodies overnight at 4°C with gentle rocking. Blots were washed with TBST, then incubated with secondary antibodies diluted in TBST + 5% BSA for 1 hr at RT with gentle rocking. Membranes were washed again with TBST, developed using SuperSignal West Femto Maximum Sensitivity Substrate (Thermo), and imaged on a ChemiDoc MP imaging system (Bio-Rad).

## Flow cytometry

After dissociating cells from culture dishes, they were pelleted and resuspended in Dulbecco’s phosphate-buffered saline (Gibco 14190-144) supplemented with 5% FBS. Forward scatter (FSC), side scatter (SSC), BFP fluorescence (BV421), YFP fluorescence (FITC), and mCherry fluorescence (PE.Texas.Red) were measured for 10,000 cells in each sample using a BD FACS Symphony or Fortessa instrument.

## Computational analyses

Pre-processing steps for high-throughput sequencing analyses were implemented as Snakemake<sup>[60](#ref-Koster2012)</sup> workflows run within Singularity containers on an HPC cluster. `Python` (v3.9.15) and `R` (v4.2.2) programming languages were used for all analyses unless mentioned otherwise. Analysis of RNA-Seq, ribosome profiling, and Perturb-Seq data, as well as gene ontology analyses are described in Supplementary Methods.

## Barcode to insert assignment

Raw data from insert-barcode linkage sequencing are in `FASTQ` format. Barcode and sgRNA insert sequences were extracted from corresponding reads and counted using `awk`; sgRNA inserts and corresponding barcodes were omitted if the sequenced sgRNA insert was not present in the designed sgRNA library (oAS1899). The remaining barcodes were aligned against themselves by first building an index with `bowtie2-build` with default options and then aligning using `bowtie2` with options `-L 19 -N 1 --all --norc --no-unal -f`. Self-alignment was used to exclude barcodes that are linked to distinct inserts or ones that are linked to the same insert but are aligned against each other by `bowtie2` (presumably due to sequencing errors). In the latter case, the barcode with the lower count is discarded in `filter_barcodes.ipynb`. The final list of insert-barcode pairs with a minimum of 5 reads is written as a comma-delimited `.csv` file for aligning barcodes from genomic DNA and mRNA sequencing below.

## Barcode counting in genomic DNA and mRNA

Raw data from sequencing barcodes in genomic DNA and mRNA are in `FASTQ` format. Barcode and UMI sequences were extracted from corresponding reads, counted using `awk`, and assigned to reporters based on their unique 6xN identifier. Only distinct barcode-UMI combinations where the barcode is present in the filtered barcodes `.csv` file from linkage sequencing are retained. The number of UMIs per barcode and associated insert are written to a `.csv` file for subsequent analyses in R. Only barcodes with a minimum of 20 UMIs were used for analysis. Barcode counts from pairs of samples were used to run MAGeCK<sup>[25](#ref-Li2014)</sup> with `--additional-rra-parameters` set to `--min-number-goodsgrna 3`. sgRNAs without a minimum of 20 UMI in one of the compared samples were set to 20 UMI counts before running MAGeCK.

## Chemicals

ISRIB (SML0843) was obtained from Sigma. Homoharringtonine (FH15974) was sourced from Biosynth. Anisomycin (A50100), hygromycin B (H75020), and puromycin dihydrochloride (P33020) were purchased from Research Products International. Nilotinib (A8232) was acquired from Apex Bio. Vemurafenib (S1267) was obtained from Selleckchem. SMG1i (HY-124719), GCN2iB (HY-112654), and eFT226 (HY-112163) were all purchased from Medchemexpress. Finally, 4E1RCat (S7370) was obtained from Selleckchem.

## Antibodies

The p38 MAPK antibody (8690, RRID: [AB_10999090](https://www.antibodyregistry.org/AB_10999090)) was obtained from Cell Signaling Technology. The phospho-p38 (Thr180/Tyr182) antibody (690201, RRID: [AB_2801132](https://www.antibodyregistry.org/AB_2801132)) was sourced from BioLegend. The GCN1 antibody (A301843AT, RRID: [AB_1264319](https://www.antibodyregistry.org/AB_1264319)) was purchased from Bethyl. The *β*-actin antibody (A5441, RRID: [AB_476744](https://www.antibodyregistry.org/AB_476744)) was obtained from Sigma. Secondary antibodies included goat anti-rabbit IgG (H+L)-HRP conjugate (1721019, RRID: [AB_11125143](https://www.antibodyregistry.org/AB_11125143)) and goat anti-mouse IgG (H+L)-HRP conjugate (1721011, RRID: [AB_11125936](https://www.antibodyregistry.org/AB_11125936)), both sourced from Bio-Rad.

## Statistics and Reproducibility

All protein immunoblotting and DNA electrophoresis experiments were repeated at least 3 times and representative gel images are shown.

# Data Availability

All high throughput sequencing data are publicly available in the NCBI SRA database under BioProject PRJNA1059490. SRA accession numbers with sample annotations are provided as Supplementary Table 5. All other data are publicly available at <https://github.com/rasilab/nugent_2024>.

# Code Availability

All software used in this study are publicly available as Docker images at <https://github.com/orgs/rasilab/packages>. Data analysis and visualization code are publicly available at <https://github.com/rasilab/nugent_2024>. Information not included in the manuscript can be publicly requested at <https://github.com/rasilab/nugent_2024/issues>.
