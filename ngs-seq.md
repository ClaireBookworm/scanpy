# Next Generation Sequencing

back then, millions of clones in 9 months for billions of dollars -> human genome project 

1. start with millions of copies (of whatever format)
2. break it in to fragments and sequence each
3. sequencer (no locational data) -> turned into "reads", the first 35-500 bps 


Platforms:
- illumina (dominates) 
- 454
- roche
- applied biosystems

**general next gen sequencing** Relies on the construction of libraries of DNA fragments that represent the DNA of the genome - instead of using bacterial cells to generate these libraries, they are made using PCR amplification of billions of DNA fragments, each attacked to a solid support. 
- pcr generated copies remain bound in proximity to the original DNA fragment (where each cluster contains about 1k identical copies of a small bit of the genome) 
- clusters sequenced in parallel 


## Illumina/Solexa
- based on dideoxy method 
- nucleotide attached to removable fluorescent molecule (diff color for each base) and a specific *chain-terminating* chemical adduct (instead of 3'-oh they have a group that blocks elongation by DNA polymerase) 
- the 4 fluorescently labeled nucleotides are added to DNA clusters on a slide, only the appropriate nucleotide (complementary to template) is covalently incorporated and camera finds the color flash
- billios of sequencing reactions are carried our simultaneously 

Assembling a whole genome you'd need the whole lane -- a basic unit you can easily separate out so you know where I came from. 

- 8 lanes
- ~160m short reads (50-70 bp) per lane (a long file with many short sequences and reads)

Use the DNA as barcodes -> knows what samples it came from

**4 basic tests**
1. sample prep
2. cluster generation
3. sequencing 
4. data analysis

### Sample Prep
add adapters to the ends of dna fragments and through reduced cycle amplification, additional motifs are introduced like indices, sequencing binding sites, and complementary regions.

### Clustering
Each fragment is isothermally amplified - flow cell (glass slide with lane, each lane has a lawn made of 2 types of oligos), where hybridization happens with one of the oligos (complementary to the adaptor on one of the fragment) and polymerase makes a strand that is kept, and the original template washed away

Strand folds over and binds to the other oligo in the flow cell, and polymerases make a double stranded bridge, and the two strandes separate, making two strands tethered to the flow cell. *Clonal amplification*, and only the forward strands are kept (3' ends blocked)

### Seqeuncing
Extension of the sequencing primer where you add nucleotides and you measure the fluroescence. 

Sequencing by sequencing - # of cycles determines the length of the read. Two reads (iintial and index read) are formed and then the 3' end is deprotected, so index 2 can be read similarly too to form double stranded bridge. 

We essneitally double read a double strand (initial 1, index 1, index 2, inittial 2)

### Data Analysis
Reads with similar fluorescnes are clustered and combined, to form a contiguous sequence.

Genomic data is easily shared and used through simple files.

## Ion Torrent Sequencing 
genome fragmented and individual fragments are attached to microscopic beads