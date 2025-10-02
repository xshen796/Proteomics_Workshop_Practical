MQ-DATAMIND ECR workshop: Proteomics and Mental Health (tutorial - local
version)
================
X Shen
02 October, 2025

Welcome! This is the repository for the MQ-DATAMIND workshop,
‘Proteomics and Mental Health’ (Practical Session II)

# Aim of the session

We will use two-sample Mendelian randomisation to test causal
relationship between protein abundance and bipolar disorder. Through
this tutorial, you will gain hands-on experience in handling and
performing quality check on pQTL and other types of GWAS sumstats. You
will also learn how to run and interpret Mendelian Randomisation
findings, and how to consolidate causal evidence using colocalisation
analysis.

We will focus on one protein, FOLH1 (Folate hydrolase 1, also known as
prostate-specific membrane antigen, PSMA), which has been implicated in
bipolar disorder by a recent Mendelian randomisation study published in
JAMA Psychiatry (Bhattacharyya U. et al. Circulating Blood-Based
Proteins in Psychopathology and Cognition: A Mendelian Randomization
Study). We will use FOLH1 as an example to illustrate the analysis
workflow, followed by functional lookup using OpenTargets. You can then
apply the same workflow to other proteins of your interest.

# Data used in this tutorial and relevant readings

------------------------------------------------------------------------

We will use publicly available summary statistics from the [Psychiatric
Genomic
Consortium](https://pgc.unc.edu/for-researchers/download-results/) and
pQTL data from [deCODE](https://www.decode.com/summarydata/)
(Ferkingstad, E. et al. Large-scale integration of the plasma proteome
with genetics and disease). Note that the complete pQTL data requires a
simple access sign up. You can download your own version for your
analysis following the link, or you could follow this tutorial using a
local copy of pre-processed pQTL data.

Readings:

-   **pQTL**: Eldjarn, G.H., Ferkingstad, E., Lund, S.H. et
    al. Large-scale plasma proteomics comparisons through genetics and
    disease associations. Nature 622, 348–358 (2023).
    <https://doi.org/10.1038/s41586-023-06563-x>

-   **Bipolar GWAS**: O’Connell, K.S., Koromina, M., van der Veen, T. et
    al. Genomics yields biological and phenotypic insights into bipolar
    disorder. Nature 639, 968–975 (2025).
    <https://doi.org/10.1038/s41586-024-08468-9>

-   **Mendelian Randomisation analysis**: Bhattacharyya U, John J, Lam
    M, Fisher J, Sun B, Baird D, Burgess S, Chen CY, Lencz T.
    Circulating Blood-Based Proteins in Psychopathology and Cognition: A
    Mendelian Randomization Study. JAMA Psychiatry. 2025 May
    1;82(5):481-491. <https://doi.org/10.1001/jamapsychiatry.2025.0033>.

# Data preparation

------------------------------------------------------------------------

## Clone this github repo

Use the script below to clone this github repository on your local
machine:

``` bash
git clone https://github.com/xshen796/Proteomics_Workshop_Practical.git
```

## Set up your R environment to access OpenGWAS data

💡See a quick guide here if you haven’t done so:
[URL](https://github.com/xshen796/Proteomics_Workshop_Practical/blob/main/Session_ii/Setup_APItoken.md)

## Data preparation

**See
[here](https://github.com/xshen796/Proteomics_Workshop_Practical/blob/main/prep/data_prep.md)
if you need more information on how we accessed and prepared the data
for this tutorial.**

# Two-sample MR: finding causal proteins to schizophrenia using cis pQTL

------------------------------------------------------------------------

## Load R packages

``` r
library(rmarkdown)
library(dplyr)
library(data.table)
library(pbapply)
library(readr)
library(Hmisc)
library(here)
library(knitr)
library(kableExtra)
library(TwoSampleMR)
library(coloc)
library("locuszoomr")
library(EnsDb.Hsapiens.v75)
```

## Load data

Let’s first load bipolar disorder GWAS sumstats and pQTL data

``` r
dat.bp_gwas = read_tsv('https://storage.googleapis.com/mhp-proteomic-sumstats/bp_2024.tsv.gz')

dat.pqtl_folh1 = read_tsv('https://storage.googleapis.com/mhp-proteomic-sumstats/Proteomics_SMP_PC0_5478_50_FOLH1_PSMA_10032022_lo_annot.txt.gz')
```

Let’s load protein annotation and see which proteins we will be looking
at

``` r
ref.somalogic = read_tsv(here::here('utils/ref_somalogic.tsv')) %>%
    dplyr::select(SeqId,UniProt,Protein_full_name,Gene,Ensembl_Gene_ID)

knitr::kable(ref.somalogic, "pipe")
```

| SeqId    | UniProt | Protein_full_name                    | Gene | Ensembl_Gene_ID |
|:---------|:--------|:-------------------------------------|:-----|:----------------|
| 3045_72  | P21246  | Pleiotrophin                         | PTN  | ENSG00000105894 |
| 9950_229 | P18627  | Lymphocyte activation gene 3 protein | LAG3 | ENSG00000089692 |

Here in this table:

-   SeqID: unique protein identifier of SOMALogic assay

-   UniProt, Protein_full_name: Uniprot ID and full name of the protein

-   Gene, Ensembl_Gene_ID: Gene symbol and Ensembl ID

## QC and prepare pQTL data (protein as exposure)

➡️ We need to select the cis region of ***FOLH1***. Find genomic
locations for ***FOLH1***
[here](https://www.genecards.org/cgi-bin/carddisp.pl?gene=FOLH1&keywords=FOLH1)

Our genome build is **GRCh37/hg19**.

![](FOLH1.jpg)

``` r
gene_boundary = data.frame(lower=49166644,upper=49230154)
gene_chr = 11
```

Select cis region sumstats: 1Mb extended area of the gene region

``` r
dat.cis_pqtl_folh1 = dat.pqtl_folh1 %>% 
  dplyr::filter(CHR==gene_chr,BP > (gene_boundary$lower - 1000000), BP < (gene_boundary$upper + 1000000))
```

Double check if data only includes cis region

``` r
range(dat.cis_pqtl_folh1$BP)
```

    ## [1] 48167081 50230092

➡️ Keep common variants.

``` r
dat.cis_pqtl_folh1 = dat.cis_pqtl_folh1 %>% dplyr::filter(eaf < 0.995, eaf > 0.005)
```

➡️ Keep variants that also present in the bipolar disorder GWAS
sumstats.

``` r
dat.cis_pqtl_folh1 = dat.cis_pqtl_folh1 %>% dplyr::filter(SNP %in% dat.bp_gwas$SNP)
```

➡️ Select genome-wide significant genetic instruments (pval \< 5e-8).

``` r
dat.cis_sig_pqtl_folh1 = dat.cis_pqtl_folh1 %>% 
  dplyr::filter(pval<5e-8)
```

➡️ Clump data using the 1000 Genome central European genotype data as
reference.

``` r
dat.cis_sig_clump_pqtl_folh1 = dat.cis_sig_pqtl_folh1 %>% 
  clump_data(.,clump_kb=1000,clump_r2=0.001) %>%
  as.data.frame
```

If you see any error message, it is likely that you have not set up your
R environment correctly. See this
[tutorial](https://github.com/xshen796/Proteomics_Workshop_Practical/blob/main/Session_ii/Setup_APItoken.md).

➡️ Reformat exposure data (pQTL)

``` r
dat.exposure.folh1 = format_data(dat.cis_sig_clump_pqtl_folh1,type='exposure')
```

## QC and prepare schizophrenia sumstats (schizophrenia as outcome)

``` r
dat.bp.outcome = dat.bp_gwas %>% data.frame %>%
  format_data(.,type='outcome') 
```

## Harmonise exposure and outcome data

Harmonise exposure and outcome data. Select the top SNP for Wald ratio
test.

``` r
dat.mr =  harmonise_data(dat.exposure.folh1,dat.bp.outcome) %>%
    arrange(pval.exposure) %>%
    head(n=1)
```

    ## Harmonising 5478_50_FOLH1_PSMA (4dOE1k) and SCZ (mTVFoE)

## Run two-sample Mendelian randomisation (Wald ratio and Steiger directionality tests)

➡️ Use the script below to run Wald ratio test. See result below:

``` r
MR.wr=mr_wald_ratio(b_exp = dat.mr$beta.exposure,
                   b_out = dat.mr$beta.outcome,
                   se_exp = dat.mr$se.exposure,
                   se_out = dat.mr$se.outcome)  

MR.wr 
```

    ## $b
    ## [1] -0.1884534
    ## 
    ## $se
    ## [1] 0.06679668
    ## 
    ## $pval
    ## [1] 0.00478297
    ## 
    ## $nsnp
    ## [1] 1

➡️ Use the script below to run MR Steiger directionality test. Here, we
use it as a complementary analysis to Wald ratio.

``` r
MR.dir = directionality_test(dat.mr)
MR.dir %>% knitr::kable(.)
```

| id.exposure | id.outcome | exposure           | outcome | snp_r2.exposure | snp_r2.outcome | correct_causal_direction | steiger_pval |
|:------------|:-----------|:-------------------|:--------|----------------:|---------------:|:-------------------------|-------------:|
| 4dOE1k      | mTVFoE     | 5478_50_FOLH1_PSMA | SCZ     |       0.0071234 |       5.81e-05 | TRUE                     |            0 |

# Colocalisation: consolidating causal evidence

------------------------------------------------------------------------

We will use the cis region pQTL data generated above for MR analyses for
colocalisation analysis. This time, we will use all the SNPs within the
cis region.

➡️ Get cis pqtl data and BP GWAS sumstats and reformat for data
harmonisation

``` r
dat.coloc.pqtl = dat.cis_pqtl_folh1 %>%
      dplyr::select(SNP,CHR,BP,effect_allele.prot=effect_allele,other_allele.prot=other_allele,
                eaf.prot=eaf,beta.prot=beta,se.prot=se,
                pval.prot=pval,samplesize.prot=samplesize,Phenotype.prot=Phenotype) %>%
      mutate(CHR=as.numeric(CHR))
```

``` r
dat.coloc.bp = dat.bp_gwas %>%
      dplyr::select(SNP,CHR,BP,effect_allele.bp=effect_allele,other_allele.bp=other_allele,
                eaf.bp=eaf,beta.bp=beta,se.bp=se,
                pval.bp=pval,samplesize.bp=samplesize,Phenotype.bp=Phenotype)
```

➡️ Merge pqtl and bp GWAS data

``` r
dat.coloc.merge = left_join(dat.coloc.pqtl,dat.coloc.bp,by=c('SNP','CHR','BP')) 
```

➡️ Align alleles

``` r
dat.coloc.merge_align = dat.coloc.merge %>% 
      mutate(ea.prot.aligned = ifelse(eaf.prot < 0.5, effect_allele.prot,other_allele.prot),
             nea.prot.aligned = ifelse(eaf.prot < 0.5, other_allele.prot,effect_allele.prot),
             maf.prot = ifelse(eaf.prot < 0.5, eaf.prot,1-eaf.prot),
             beta.prot.aligned = ifelse(eaf.prot < 0.5, beta.prot,-beta.prot)) %>% 
      mutate(ea.bp.aligned = ifelse(eaf.bp < 0.5, effect_allele.bp,other_allele.bp),
             nea.bp.aligned = ifelse(eaf.bp < 0.5, other_allele.bp,effect_allele.bp),
             maf.bp = ifelse(eaf.bp < 0.5, eaf.bp,1-eaf.bp),
             beta.bp.aligned = ifelse(eaf.bp < 0.5, beta.bp,-beta.bp)) %>% 
      mutate(flipped = ea.prot.aligned!=ea.bp.aligned) %>% 
      dplyr::filter(!is.na(pval.bp),maf.bp>0.001) %>% 
      .[!duplicated(.$SNP),]

sum.flipped = sum(dat.coloc.merge_align$flipped)
```

➡️ Format data for SuSie colocalisation analysis

``` r
D1 <- list(
             type = "quant", # quantitative trait
             beta = dat.coloc.merge_align$beta.prot.aligned,
             varbeta = dat.coloc.merge_align$se.prot^2, # note that this is standard error squared
             pvalues = dat.coloc.merge_align$pval.prot,
             N = dat.coloc.merge_align$samplesize.prot,
             MAF = dat.coloc.merge_align$maf.prot,
             snp = dat.coloc.merge_align$SNP,
             sdY = 1 # external information
            )

D2 <- list(
          type = "cc", # case-control trait
          beta = dat.coloc.merge_align$beta.bp.aligned,
          varbeta = dat.coloc.merge_align$se.bp^2, # note that this is standard error squared
          pvalues = dat.coloc.merge_align$pval.bp,
          N = dat.coloc.merge_align$samplesize.bp,
          s = 0.07,
          MAF = dat.coloc.merge_align$maf.bp,
          snp = dat.coloc.merge_align$SNP
          )
```

➡️ Run SuSie colocalisation analysis and check your results

``` r
res.coloc = coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
```

    ## Warning in check_dataset(d = dataset2, 2): minimum p value is: 2.1146e-05
    ## If this is what you expected, this is not a problem.
    ## If this is not as small as you expected, please check you supplied var(beta) and not sd(beta) for the varbeta argument. If that's not the explanation, please check the 02_data vignette.

    ## PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
    ##  1.29e-75  6.40e-02  4.00e-75  1.98e-01  7.38e-01 
    ## [1] "PP abf for shared variant: 73.8%"

``` r
res.coloc
```

    ## Coloc analysis of trait 1, trait 2

    ## 
    ## SNP Priors

    ##    p1    p2   p12 
    ## 1e-04 1e-04 1e-05

    ## 
    ## Hypothesis Priors

    ##          H0    H1    H2        H3     H4
    ##  -0.4542461 0.549 0.549 0.3013461 0.0549

    ## 
    ## Posterior

    ##        nsnps           H0           H1           H2           H3           H4 
    ## 5.490000e+03 1.289080e-75 6.403637e-02 4.003058e-75 1.981181e-01 7.378455e-01

➡️ Visualise the results

``` r
dat.plot.folh1 = dat.coloc.pqtl %>% dplyr::select(chrom=CHR,pos=BP,rsid=SNP,other_allele=other_allele.prot,effect_allele=effect_allele.prot,p=pval.prot,beta=beta.prot,se=se.prot)

loc <- locus(gene = 'FOLH1', dat.plot.folh1, flank = 1e5,ens_db = "EnsDb.Hsapiens.v75")

locus_plot(loc)
```

![](FOLH1_Bipolar_practical_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#detach("EnsDb.Hsapiens.v75")
```

**pQTL in FOLH1 genomic region**

``` r
dat.plot.folh1 = dat.coloc.pqtl %>% dplyr::select(chrom=CHR,pos=BP,rsid=SNP,other_allele=other_allele.prot,effect_allele=effect_allele.prot,p=pval.prot,beta=beta.prot,se=se.prot)

loc <- locus(gene = 'FOLH1', dat.plot.folh1, flank = 1e6,ens_db = "EnsDb.Hsapiens.v75")

locus_plot(loc)
```

![](FOLH1_Bipolar_practical_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

**Bipolar disorder QTL in FOLH1 genomic region**

``` r
dat.plot.bp = dat.coloc.bp %>% dplyr::select(chrom=CHR,pos=BP,rsid=SNP,other_allele=other_allele.bp,effect_allele=effect_allele.bp,p=pval.bp,beta=beta.bp,se=se.bp)

loc <- locus(gene = 'FOLH1', dat.plot.bp, flank = 1e6,ens_db = "EnsDb.Hsapiens.v75")

locus_plot(loc)
```

![](FOLH1_Bipolar_practical_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## Interpreting the findings

### Gene lookup (GeneCards)

➡️ Go to [GeneCards](https://www.genecards.org/)

### Mapping FOLH1 to drugs

-   Source 1: OpenTargets

-   Source 2: DrugBank online
