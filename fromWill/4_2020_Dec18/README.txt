This directory contains files for the secondary structure and sequence evolution analysis of PNMA (pnma.subtree/)and Mart (mart.subtree) family genes.

Columns are comma-delimited, and files are delimited with > (for ease of parsing) 

Notes: Alignments could be trimmed/some sequences deleted for better structural analysis - would this help or hinder phylogenetic analysis? 

Next steps: Trim these alignments to contain only the most interesting sequences (from structural/evolutionary perspective), re-align. 

file, num.seqs,	description
>2020.12.14.kingdom.vertebrates.canonical.retrotransposon.mart.pnma.subtree.png,0, A schematic of the tree and phylogenetic neighbors of PNMA and Mart family genes 
>mart.subtree/2020.12.10.CA.full.seqs.mart.neighborhood.extended.human.only.faa, 24, Full-length human sequences with HMM hits for CA domain (including those excluded by downstream phylogenetic analysis)(8 included downstream - RTL1|RTL5|6 isoforms of PEG10)(16 excluded from downstream) 
>mart.subtree/2020.12.10.CA.only.seqs.refseq.mart.extended.faa, 24, CA domain only of human sequences with HMM hits for CA domain (including those excluded by downstream phylogenetic analysis)(8 included downstream - RTL1|RTL5|6 isoforms of PEG10)(16 excluded from downstream)-Not used for alignments to allow for full-length sequence observations
>mart.subtree/2020.12.11.CA.full.mart.neighborhood.extended.Ty3.faa, 119, Full length sequences from the Mart family subtree and excluded human sequences (extended) (upper in schematic, and HHpred-predicted closest capsid of Ty3 (n=2 Ty3 sequences - Uniprot and Repbase)(n=32 human sequences - some duplicates because of concatenation strategy used to combine subtree family sequences and human sequences - these human sequences are reduced to n=14 by cdhit)(n=85 vertebrate Gypsy retrotransposons)(NOTE: prefix for these sequences reduced to 2020.12.*.CA.full.mart.neighborhood.* in downstream workflow)
>mart.subtree/2020.12.11.CA.full.mart.neighborhood.cdhit.t1.linsi.aln, 101, MAFFT-alignment of extended Mart neighborhood (L-INS-I strategy)(Note: CD-Hit removed 18 identical sequences, all human sequences prior to alignment)(NOTE: This file was used to build espript file)
>mart.subtree/2020.12.11.CA.full.mart.neighborhood.cdhit.t1.linsi.trimal.gt0.8.aln, 101, Trimmed MSA (all columns with sequence in less than 80% of the rows removed)(NOTE: this file was used for Ali2d secondary structure prediction)
>mart.subtree/2020.12.13.CA.full.mart.neighborhood.cdhit.t1.linsi.espript.pdf, 101, Espript-prepared figure of alignment across the helices of the CA domain with Ty3 structures overlaid the alignment (NTD top - CTD bottom)
>mart.subtree/2020.12.13.CA.full.mart.neighborhood.cdhit.t1.linsi.trimal.gt0.8.ali2d.out, 100, Ali2D output (I think this can be easily manipulated for a figure in R, but I haven't played with it yet)(NOTE: Repbase Ty3 sequence removed because Ali2D has a 100-sequence limit) 
>mart.subtree/2020.12.13.CA.full.mart.neighborhood.cdhit.t1.linsi.trimal.gt0.8.ali2d.out.01.png, 46, Top half of screenshot of Ali2d visual output - shows conserved helices across Gypsy neighbors and human sequences. Also shows that some human genes are missing conserved regions. (NOTE: several residues of conserved zinc-finger (CCHC motif) at the far right (CXXC is only one visible) (NOTE: RTL1 and RTL5 are the first and third sequences, and lack this motif. RTL3 -  NP_689907.1 - appears to have the motif, but was presumably excluded in the phylogenetic analysis because of several gappy regions in the alignment. I think it warrants more exploration) 
>mart.subtree/2020.12.13.CA.full.mart.neighborhood.cdhit.t1.linsi.trimal.gt0.8.ali2d.out.02.png, 48, Second half of Ali2d visual output. Note: Sequences 47-50 are not visible, but the patterns seen in other Gypsys hold). 
>pnma.subtree/2020.12.14.CA.full.seqs.refseq.pnma.extended.faa, 32,  Full length human sequences with HMM hits for CA domain (including those excluded by downstream phylogenetic analysis)(25 included in downstream analysis)(7 excluded from downstream analysis)
>pnma.subtree/2020.12.14.CA.only.seqs.refseq.pnma.extended.faa, 32, CA domain only of human sequences with HMM hits for CA domain (including those excluded by downstream phylogenetic analysis)
>pnma.subtree/2020.12.15.CA.full.seqs.neighborhood.extended.dARCs.faa, 79, Full length sequences from the PNMA family subtree and excluded human sequences (extended) (lower tree in schematic, and HHpred-predicted closest capsid of dARC (n=2 dARC sequences - Uniprot dARC1 and Uniprot dARC2)(n=55 human sequences - some duplicates because of concatenation strategy used to combine subtree family sequences and human sequences - these human sequences are reduced to n=17 by cdhit)(n=23 vertebrate Gypsy retrotransposons)
>pnma.subtree/2020.12.15.CA.full.seqs.neighborhood.extended.dARCs.cdhit.t1.linsi.aln, 42, MAFFT-alignment of extended PNMA neighborhood (L-INS-I strategy)(Note: CD-Hit removed 38 identical sequences, all human sequences prior to alignment)(NOTE: This file was used to build espript file)
>pnma.subtree/2020.12.15.CA.full.seqs.neighborhood.extended.dARCs.cdhit.t1.linsi.espript.pdf, 42, Espript-prepared figure of alignment across the helices of the CA domain with dARC structures overlaid the alignment (dARC1 top - dARC2 bottom)
>pnma.subtree/2020.12.15.CA.full.seqs.neighborhood.extended.dARCs.cdhit.t1.linsi.trimal.gt0.8.aln, 42, Trimmed MSA (all columns with sequence in less than 80% of the rows removed)(NOTE: this file was used for Ali2d secondary structure prediction)
>pnma.subtree/2020.12.15.CA.full.seqs.neighborhood.extended.dARCs.cdhit.t1.linsi.trimal.gt0.8.ali2d.out, 42, Ali2D output (I think this can be easily manipulated for a figure in R, but I haven't played with it yet)
>pnma.subtree/2020.12.15.CA.full.seqs.neighborhood.extended.dARCs.cdhit.t1.linsi.trimal.gt0.8.ali2d.png, 42, Ali2d visual output - shows conserved helices across Gypsy neighbors and human sequences. Also shows that some human genes are missing conserved regions. (NOTE: the CCHC motif is not visible in the far right, it is trimmed - however it is visible in espript. I think it warrants more exploration, and a next step is to delete all sequences that lack it and re-align) 











 

