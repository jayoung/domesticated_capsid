library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(org.Hs.eg.db)

# library(BSgenome.Hsapiens.UCSC.hg38)


# library(phastCons100way.UCSC.hg38) # don't use phastcons - they are the smoothed, not the basewise scores




PNMA3_region <- GRangesForUCSCGenome("hg38", "chrX", IRanges(start=153055486, end=153061457))

 
tracks["Conservation"]
# "cons100way" 

query <- ucscTableQuery(session, "Conservation", PNMA3_region)

## list the table names
tableNames(query)
# [1] "phyloP100way"            "phastCons100way"        
# [3] "phastConsElements100way" "multiz100way"           


## get the phastCons30way track
tableName(query) <- "phyloP100way"
track(query)  # a GRanges object, each region 1bp wide

## vignette shows using rtracklayer to get multiz track, but it doesn't actually get mafs




# https://www.bioconductor.org/help/course-materials/2015/UseBioconductorFeb2015/A01.5_Annotation.html



TxDb.Hsapiens.UCSC.hg38.knownGene
# UCSC Track: GENCODE v32
# Type of Gene ID: Entrez Gene ID
### this GENCODE track has a LOT of isoforms


session <- browserSession("UCSC")
genome(session) <- "hg38"

tracks <- trackNames(session)


## want cds by transcript so as to be able to pick out isoforms
hg38_exons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx")
hg38_cds <- cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx")
    # names of those are internal IDs for the package, but I can get the ENST names using id2name

###
# hg38_exons2 <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx", use.names=TRUE)
# hg38_cds2 <- cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "tx", use.names=TRUE)
# Warning message:
# In .set_group_names(grl, use.names, txdb, by) :
#   some group names are NAs or duplicated


txIDs <- id2name(TxDb.Hsapiens.UCSC.hg38.knownGene, feature.type="tx")
length(txIDs)
# [1] 247541
length(unique(txIDs))
# [1] 247380

txIDs[ which(duplicated(txIDs)) ]

# this is on the X and Y chromosome
txIDs[ which(txIDs=="ENST00000381223.9") ]
#              223078              227007 
# "ENST00000381223.9" "ENST00000381223.9" 

hg38_cds[["223078"]]
hg38_cds[["227007"]]


## example gene:  NCBI gene ID = 1 = A1BG
## example gene:  NCBI gene ID = 10 = NAT2.  2 isoforms, each with 2 exons

## example gene:  NCBI gene ID = 60598 = KCNK15. 1 isoform, 2 coding exons
hg38_exons[["60598"]]
hg38_cds[["60598"]]
    # it just shows the two exons

## example gene: 
# PNMA3  NCBI geneID=29944 chrX:153,055,486-153,061,457
# two Refseq isoforms, encoding different proteins: NM_013364.6, NM_001282535.2 
# ensembl gene: ENSG00000183837.9
# three ensembl transcripts (one subject to NMD):
#     PNMA3-203 = ENST00000619635.1 = 3627bp 455aa, Protein coding
#     PNMA3-202 = ENST00000593810.2 = 1432bp 463aa, Protein coding
#     PNMA3-201 = ENST00000424805.1 = 3158bp 463a, Nonsense mediated decay (not shown on UCSC website- it shows the 'basic set' rather than the 'comprehensive' set)


hg38_exons[["29944"]]
hg38_cds[["29944"]]
# 

### xxx what transcripts go with gene PNMA3 / 29944

id2name(TxDb.Hsapiens.UCSC.hg38.knownGene, feature.type="cds")[1:40]


# get entrez gene ID from gene name
select(org.Hs.eg.db, "PNMA3", "ENTREZID", "SYMBOL")[["ENTREZID"]]
# [1] "29944"

# get transcript name(s) from Entrez gene identifier

select(TxDb.Hsapiens.UCSC.hg38.knownGene, "29944", "TXNAME", "GENEID")[["TXNAME"]]
# [1] "ENST00000619635.1" "ENST00000424805.1"
## why is ENST00000593810.2 not shown? Ensembl lists it as transcript support level not analyzed - that might be why.


select(TxDb.Hsapiens.UCSC.hg38.knownGene, keys="29944", columns=c("GENEID", "TXNAME", "TXID"), keytype="GENEID")
#   GENEID   TXID            TXNAME
# 1  29944 222719 ENST00000619635.1
# 2  29944 222720 ENST00000424805.1

select(TxDb.Hsapiens.UCSC.hg38.knownGene, keys="ENST00000593810.2", columns=c("GENEID","TXNAME"), keytype="TXNAME")
# 'select()' returned 1:1 mapping between keys and columns
#             TXNAME GENEID
# 1 ENST00000593810.2   <NA>

txIDs[ which(txIDs %in% c("ENST00000619635.1", "ENST00000424805.1")) ]
#              222719              222720 
# "ENST00000619635.1" "ENST00000424805.1" 

txIDs[ which(txIDs %in% c("ENST00000593810.2")) ]
#              222721 
# "ENST00000593810.2" 

## these are the exons of the CDS of the three transcripts for PNMA3
hg38_cds[["222719"]]
hg38_cds[["222720"]]
hg38_cds[["222721"]]
    # the last two are the same as each other

### xxx how would I get the conservation tracks for the exons in one of those transcripts? (ENST00000619635.1):

hg38_cds[["222719"]]

