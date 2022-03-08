# Single Cell Sequencing notes (USCF)

- single cell sequencing (RNA)
- plate based methods 
- microfluidic-based methods
- combinatorial indexing methods
- other types of single cell sequencing methods (non-RNA)


## Bulk vs. single cell analogy 
blood has a bunch of different forms of cells -> all the different tyeps means that having bulk means you can't find the characteristics of that specific cell (you mix it up instead)

**Single cell data can identify cell subtypes**

(looking at the graphs from `paga-paul15` you can remember how there were different labels and colors)
-> each dot represents a single cell, and clustered means the closer they are the most similar they are 

(*zhang et al 2016* interesting paper to read?)
SCS can be further broken dwown more and more as opposed to flow cytometry

scRNA-seq output has increased significantly

## plate-based SMART-seq
uses the oligo-dt (poly-A tail) and reverse transcriptase binds to that to transcript -> and then the RT adds some unpaired C bases that binds to some G bases -- oligo combination

ends up with cDNA with **pcr handles** on both sides -- happens multiples times 
*parallelized* the process! 

## DropSeq 
uses microfluidics (different methods of creating very small volume droplets and doing mole bio manipulations in them)

- chip with beads -> cells flow in and mixed together -> oil flows in and causes aqueous phases -> finally very small amt of droplets will have both droplet and bead  (bead has different barode so you can capture them)
- RT reaction with the template switching process (like in SMART) and all the sample cells are stuck to one bead 
- PCR 
- sequencing and analysis (mRNA is mapped to its cell of origin and gene of origin, each cells pool of mRNA can be analyzed)

this doesn't require an insane amt of wells 

**POISSON** DISTRBUTION

issue with this method is that most of the droplets don't have anything and only a fraction of droploets even have a cell

you can increase the [] of beads and cells, but you would have a doublet rate that's very high 

> doublets are bad because tyou have multiple cells per droplet/bead

## 10x Genomics 

barcoded beads + cell/reagents -> oil -> *emulsion*

- most of the droplets have beads so >90% have one bead and only one bead 
- in drops there are already RT within the droplets 
- at the end you break the emulson, amplify cDNA, and construct the library to get the sequence 

only takes 15 minutes and the chip can do 200k cells at once (breaking emultion takes half a day)

## Microwells - parallel scs without droplets 
- seqwell
- celsee genesis
- becton dickinson rhapdosy

there are very small microwell array (hundreds of thousands of wells in a chip) and you laydown the barcoded beads (one and only 1 bead in a well) and you apply cells on top  (blue colored bead)

cells settle down in wells -> lysis and enzymatic reactions -> attack the cells to get the RNA and complete RT 

## Combinatorial Indexing 
- single cell conbinatorial indexing and SPLIT-seq
- *in situ* reactions that add barcodes 
- split pooling in between each step 

**Fixed** cells -> sort a few dozen cells into each well -> RT with barcoded oligo-DT (each well has different sequence) -> combine thecells and shuffle them up and redistribute them into a new set of wells 

Taking the cDNA in each of the wells that adds a *second* barcode -- chances are that each cell came from a different well in the first place . This has 10k combinations and the finally library is:

`RT | UMI | olido(dT) | cDNA`

- you can have multple plates that are larger or do it multople multiple times, so you have millions of different combinations
- almost impossible for them to be in the same well again and also means you can sequence so so many cells at the same time.

## Quantifying proteins - CITE-seq
antibodies bind to proteins of interest (instead of fluorophore) they bind to a DNA barcode for each andibody + the polyA tail (which allows capture by the SCS methods)

This helps tell us how many **copies** of this antibody was found in the cell. You can also figure out what genes are turned on and off (not entirely sure how?)

## Multiplets limit laoding concentration
Doublets are bad (want to limit) and so are multiplets. There are many methods to reduce the chances of it being in the data at the end.

###  demuxlet - Multiplexing by genotype 
- encapsulate geenticallly diverse cells in drops 
- cDNA labeled with genetic and cell barcode 

Looking at each individual cell, we see a specific genotype -> can then find the donor and what patient it's from.

If you have a *doublet*, it's likely that the two cells came from two patients/and have different doublets, so if we see that we remove from the dataset (has to be genetically distinct so might not work for everything).

### Non-genetic multiplexing strategies
- samples -> label wtih hashtag oligos (HTO)
- cell pooling 

pools the cells after labeling the antibodies with barcodes -- so with 2+ cells you can detect bc there are 2 different barcodes 

- lipidated DNA that embeds in the DNA of the cell
- depends on accessible membrane in the cell 

## Conclusions 
- many types of RNA, DNA, protein
- generations high-dimensional data (20k genes) -> reveals biology hidden from lower dimensional techniques 
- lower thruput and more *expensive* compared to flow / mass cytometry
- microscopy methods (doesn't need sequencing)
	- seqFISH, MERFISH, CODEX
	- lower thruput but provide spatial information


# Science SCS Introduction

