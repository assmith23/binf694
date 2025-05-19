#
## A Manning Smith
## 5/20/2025



# Introduction

A single nucleotide variant (SNP), a single nucleotide subsitution at a specific genome postion, are the most common sequence variations [**[#]**](#sherry1999).

# Scoreing Methods
## PolyPhen
The PolyPhen-2 is a tool utilized for the annotation of SNPs. The foundation of the classifer is a "machine-learning method" [**[#]**](http://genetics.bwh.harvard.edu/wiki/!pph2/about). The tool works to define various amino acid replacements, by feeding "structure-based features of the the subsitution site to a probabilistic classifer" [**[#]**](http://genetics.bwh.harvard.edu/wiki/!pph2/overview).

As a probability value that range of the score is from $0-1$ being that $0$ has less of a probability of being damaging vs values closer to $1$ have a higher probability of being damaging. The performace of the tool is said to achieve a true postivite prediction of "92% and 73% on HumDiv and HumVar, respectively" [**[#]**](#Adzhubei2010).

## GranthamScore
The Grantham score is a prediction, in an evolutionary sense,  of the distance between two amino acids [**[#]**](https://ionreporter.thermofisher.com/ionreporter/help/GUID-D9DFB21C-652D-4F95-8132-A0C442F65399.html). The score takes into account the composition, polarity, and molecular volume [**[#]**](#grantham1974). This score is a more straight forward calculation as seen below.

$$
    D_{ij} = \Big[\alpha (c_i-c_j)^2 + \Beta (p_i-p_j)^2 + \rou(v_i-v_j)^2 \Big]^{1/2}
$$

The score ranges from $5$ to $215$. A lower score, closer to $5$, indicates a less significant subsitution, while a high score represents that of a more critical subsitution. It uses a matrix, $D$, to related the distances between amino acids, thus representing the severity of change between two amino acids.

## CADD
Combined Annotation Dependent Depletion (CADD) is a tool used to determine the harmfulness of a SNP. Different from our other two scores, CADD takes into account variants throughout the entire genome. The values is calculated via a machine learning model, "trained on a binary distinction between simulated de novo variants and variants that have arisen and become fixed in human populations..." [**[#]**](#philipp2019). The value combines aspects like, conservation metrics, transcription factor binding, and genomic and epigenomic annotations, along with a few other parameters.

The score is commonly defined as "scaled C-scores" described as the "rank of each variant relative to all possible 8.6 billion substitutions" [**[#]**](https://cadd.gs.washington.edu/info). It ranges from $10-99$, represented as $10 * -\log{\text{rank}}$. A score above 20, is "predicted to be amoung the $1.0%$ most deleterious possible subsitutions" in the genome [**[#]**](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&chr=chrX&g=caddSuper). 

## Score Analysis

[**Table 1**](#tab1) highlights the summary statistics of our variant observations for the three scores highlighted. There were a large amount of observations removed for having missing values. CADD had significantly less amount of observations removed, likly due to the calculation being applied to non-coding regions, unlike the other two measures.

Our scales align with that outlined above for each score respectivly. The means of our scores represent that of less severe variants. This aligns with the health profile of the individual as the majority of the that variants are not signficant.

[**Figure 1**](#fig1) outline lines the distribution of our `complete_cases` dataset, which is comprised off the variants that have values for all three of our scores. Our scores present significantly right skewed distributions aligning with the idea that the majority of the variants are not significant.

[**Figure 2**](#fig2) outlines the correlations between the three different scores. There is a reletivtly high correlation between `PolyPhen vs CADD`, suggesting that these scores seem to agree with each other in the signifance of the variant. Our little to non correlation between `PolyPhen vs Grantham` suggest that the amino acid subsituations generally don't have a high probability of imapce. Lastly, `Grantham vs CADD` has no correlation suggesting that amino acid changes have no relationship with deleteriousness.

Based on these scoring patterns, we can now identify specific variants that highlight the need for further investigation, focusing on outliers that demonstrate potentially significant biological impact according to our methods.

# Variant Selection



[**Table 2**](#tab2) lists some basic information about the selected variants used in our analysis. 

## GALNTL5 - rsID: 6960270
The Polypeptide GalNAc transferase 15 (GALNTL5) plays an important role in the male reproductive system as it is mostly expressed in the testis and plays a role in sperm development.

### Gene Analysis

The first variant is on chromosome 7 at position 151982987. The variant is a missense variant being heterozygous with a change from T to C, changing the amino acid from CYS to ARG, at position 177/251 near the splice site.



### Molecular Impact
### Connection to Disease or Trait
### Assessment of Genotype / Phenotype

### Population Distribution

## Variant 2
### Molecular Impact
### Connection to Disease or Trait
### Assessment of Genotype / Phenotype
### Population Distribution

## Variant 3
### Molecular Impact
### Connection to Disease or Trait
### Assessment of Genotype / Phenotype
### Population Distribution

# Discussion

# Resources

<div id="Adzhubei2010"></div>

**[#]**. Adzhubei, I. A., Schmidt, S., Peshkin, L., Ramensky, V. E., Gerasimova, A., Bork, P., Kondrashov, A. S., & Sunyaev, S. R. (2010). A method and server for predicting damaging missense mutations. Nature methods, 7(4), 248–249. https://doi.org/10.1038/nmeth0410-248


<div id="philipp2019"></div>

**[#]**. Philipp Rentzsch, Daniela Witten, Gregory M Cooper, Jay Shendure, Martin Kircher, CADD: predicting the deleteriousness of variants throughout the human genome, Nucleic Acids Research, Volume 47, Issue D1, 08 January 2019, Pages D886–D894, https://doi.org/10.1093/nar/gky1016


<div id="grantham1974"></div>

**[#]**. R. Grantham ,Amino Acid Difference Formula to Help Explain Protein Evolution.Science185,862-864(1974).DOI:10.1126/science.185.4154.862


<div id="Takasaki2014"></div>

Takasaki, N., Tachibana, K., Ogasawara, S., Matsuzaki, H., Hagiuda, J., Ishikawa, H., Mochida, K., Inoue, K., Ogonuki, N., Ogura, A., Noce, T., Ito, C., Toshimori, K., & Narimatsu, H. (2014). A heterozygous mutation of GALNTL5 affects male infertility with impairment of sperm motility. Proceedings of the National Academy of Sciences of the United States of America, 111(3), 1120–1125. https://doi.org/10.1073/pnas.1310777111


<div id="sheery1999"></div>

**[#]**. Sherry, S. T., & Sirotkin, K. (1999). dbSNP—Database for Single Nucleotide Polymorphisms and Other Classes of Minor Genetic Variation. Genome Research, 9(8), 677-679. https://doi.org/10.1101/gr.9.8.677

## AI Code

<div id="ai_code1"></div>

[**[#]**](https://claude.ai/share/c65847ac-89d9-4dc0-ae30-8c8f3337cfaf). Clause was utilzied to optimzie the basic information of a list of genes to allow for further research. Clause was prompted with:
```
Research the following gene list. I want the link to the uniport entry for human. I want the gene name and the description of the gene. Save the information in a dictionary with the gene as the key. The write the python code to loop through the dictionary and save the results in an excel file.

The columns should be gen_symbol, gene_name, gene_description, uniport link.

I then want another lookup to get a excel of all the variant information on a gene. This output should be gene_symbol, ID, position, description with the rs from the dbSNP.
```
The provide code was directly used with no edits via the author. The author validated the results to ensure accuracy. There is no impact to the analysis or manucript of the paper as the results were used for further research. This code saved the user 2+ hours of data mining through UniProt.

# Tables and Figures

<div id="tab1"></div>

#### Table 1


<div id="tab2"></div>

#### Table 2


<div id="tab3"></div>

#### Table 3


<div id="fig1"></div>


#### Figure 1.0


<div id="fig2"></div>

#### Figure 2.0