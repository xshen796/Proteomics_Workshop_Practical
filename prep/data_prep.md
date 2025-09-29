pQTL, SCZ GWAS sumstats, reference data preparation
================
X Shen
29 September, 2025

# Prepare annotation file

## Files required for our analysis

-   Protein assay annotations. Key information include: UniProt ID, gene
    name, gene Ensembl IDs, SeqIDs. This is usually assay-specific.
    deCODE proteomic data used a SOMALogic assay. More information can
    be found in the
    [paper](https://www.nature.com/articles/s41588-021-00978-w). We will
    use their [Supplementary Table
    1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-021-00978-w/MediaObjects/41588_2021_978_MOESM4_ESM.xlsx).

    The protein annotation table looks this:

    | SeqId    | Protein (short name) | Protein (full name)                                | Gene   | UniProt | Organism | Type    | Ensembl.Gene.ID |
    |----------|----------------------|----------------------------------------------------|--------|---------|----------|---------|-----------------|
    | 10000_28 | CRBB2                | Beta-crystallin B2                                 | CRYBB2 | P43320  | human    | Protein | ENSG00000244752 |
    | 10001_7  | c-Raf                | RAF proto-oncogene serine/threonine-protein kinase | RAF1   | P04049  | human    | Protein | ENSG00000132155 |
    | 10003_15 | ZNF41                | Zinc finger protein 41                             | ZNF41  | P51814  | human    | Protein | ENSG00000147124 |

    A local copy is stored in utils/ref_somalogic.tsv

-   Gene annotations. Key information is start and end genomic positions
    of a given protein-encoding gene. We can query using the ‘bioMart’ R
    package. Check the tutorial
    [here](https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html).
    View online for more info
    [here](https://www.ensembl.org/info/data/biomart/how_to_use_biomart.html).

``` r
# Map to ensembl IDs
## Get reference ensembl database
mart_prot <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # human
Uniprot_ensembl = getBM(
  attributes=c('ensembl_gene_id','uniprotswissprot',"hgnc_symbol","chromosome_name","start_position","end_position"), 
  mart = mart_prot)
colnames(Uniprot_ensembl) <- c("Ensembl_ID", "UniProt","Gene_symbol","CHR","start_pos","end_pos")

knitr::kable(Uniprot_ensembl %>% arrange(CHR) %>% head(5), "pipe")
```

| Ensembl_ID      | UniProt | Gene_symbol | CHR | start_pos | end_pos |
|:----------------|:--------|:------------|:----|----------:|--------:|
| ENSG00000142611 |         | PRDM16      | 1   |   3069168 | 3438621 |
| ENSG00000142611 | Q9HAZ2  | PRDM16      | 1   |   3069168 | 3438621 |
| ENSG00000149527 |         | PLCH2       | 1   |   2425980 | 2505532 |
| ENSG00000149527 | O75038  | PLCH2       | 1   |   2425980 | 2505532 |
| ENSG00000142606 |         | MMEL1       | 1   |   2590639 | 2633016 |

``` r
write_tsv(Uniprot_ensembl,here::here('utils/uniprot_ensembl.tsv'))
```

    A local copy of Ensembl annotation is stored in utils/uniprot_ensembl.tsv

-   Variants annotation for pQTL data

    Sometimes researchers choose to store additional variant
    information, such as allele frequency, that tend to be consistent
    across multiple analyses, separately from individual GWAS sumstats.

    You can request this information for the pQTL data from
    [here](https://download.decode.is/form/2021/assocvariants.annotated.txt.gz).

    A copy of this information is stored on cloud.

# Prepare pQTL data

``` r
annot_variant = read_tsv('data/assocvariants.annotated.txt.gz') %>%
    dplyr::select(SNP=rsids,effect_allele=effectAllele,eaf=effectAlleleFreq,Name) %>%
    mutate(n_allele = nchar(effect_allele)) %>%
    filter(n_allele==1)
ref.somalogic = read_tsv(here::here('utils/ref_somalogic.tsv')) %>%
    left_join(.,read_tsv(here::here('utils/uniprot_ensembl.tsv')),by='UniProt') %>%
    mutate(CHR=as.numeric(CHR))

format_pqtl_exposure <- function(f.x.dat,x.ref){
  # Get annotation
  Symbol = f.x.dat %>% basename %>% gsub('_10032022.txt.gz','',.) %>% gsub('Proteomics_SMP_PC0_','',.) 
  x.seqid = Symbol %>% str_split(.,'_') %>% .[[1]] %>% head(2) %>% paste0(.,collapse='_')
  x.ref = x.ref %>% filter(SeqId == x.seqid)
  x.chr = x.ref$CHR
  x.start = x.ref$start_pos
  x.end = x.ref$end_pos

  # Output name
  f.out = f.x.dat %>% basename %>% gsub('.txt.gz','_lo_annot.txt.gz',.)

  # Liftover from hg38 to hg19
  x.dat.liftover = read_tsv(f.x.dat) %>%
    dplyr::filter(Chrom!='chrX') %>%
    mutate(CHR=Chrom %>% gsub('chr','',.) %>% as.numeric,
           p=10^(-minus_log10_pval)) %>%
    dplyr::select(SNP=rsids,CHR,BP=Pos,MAF=ImpMAF,A1=effectAllele,A2=otherAllele,beta=Beta,se=SE,p,N,Name) %>% 
    MungeSumstats::liftover(sumstats_dt = ., 
                            ref_genome = "hg38",
                            convert_ref_genome = "hg19",end_col='BP') %>% 
    mutate(Phenotype=Symbol) 

  # Add variant annotation
  x.dat.liftover_annot = x.dat.liftover %>%
    dplyr::select(CHR,BP,SNP,effect_allele=A1,other_allele=A2,beta,se,pval=p,samplesize=N,Phenotype,Name) %>% 
    dplyr::filter(CHR==x.chr) %>% 
    left_join(.,annot_variant,by=c('Name')) %>%
    dplyr::select(-Name)
  
  # Write preprocessed file
  dir.create(here::here('./data/pQTL_hg19_annot/'), showWarnings = FALSE)
  write_tsv(x.dat.liftover_annot,here::here(paste0('data/pQTL_hg19_annot/',f.out)))
}

# Prepare pQTL data and write files in data/pQTL_hg19_annot 

list.files(here::here('data/deCODE_pQTL/'),full.names = T) %>% 
  .[grep('.gz$',.)] %>% .[3] %>%
  as.list %>% 
  pblapply(format_pqtl_exposure,x.ref=ref.somalogic)
```

# Prepare SCZ sumstats

Download SCZ sumstats
[here](https://pgc.unc.edu/for-researchers/download-results/)

``` r
scz.gwas = read_tsv(here::here('data/daner_PGC_SCZ_w3_90_0518d_eur.gz'))
```

``` r
scz.gwas = scz.gwas %>% 
 mutate(beta = log(OR)) %>% 
 dplyr::select(SNP,effect_allele=A1,other_allele=A2,
                       eaf=FRQ_U_77258,beta,se=SE,pval=P,samplesize=Neff) %>% 
         mutate(Phenotype='SCZ') 

write_tsv(scz.gwas,here::here('data/scz_2022.tsv.gz'))
```

A preprocessed copy is stored on cloud.
