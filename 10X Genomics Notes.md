10X Genomics Notes
Starting material for scRNA-Seq is cell suspension to which all reagents from 10X's "Single Cell Gene Expression Kit" are added then run through the Chromium machine,resulting in a single-cell sequencing library.

Key steps include: Chromium microfluidic chip traps cells in oil bound droplets each of which contains a single cell and a microgel bead. This oil droplet is called a "GEM" gel bead in emulsion or GEM. 
Each microbead (of which there are more than 3 million unique ones) is coated with oligonucleotides / primers which have multiple parts to it. 

The first part is a TrueSeq read1 from Illumina needed to prepare a sequencing library. This sequence is identical for all gel beads.

The second part is a 16mer which is a 10X barcode unique to that gel bead (thus, there are more than 3 million 10X barcodes) This means, all RNA molecules from a single cell trapped in the GEM will carry this barcode. This then, identifies the cell.

the 3rd part is the UMI: a 12mer which identifies each RNA molecule individually. 12 bases randomized gives 16.7 million unique 12 mers. UMI stands for unique moelcular identifier, this allows for quantification of the RNA in each cell. Each RNA captured gets 1 amongst 16.7 million possible UMIs. So within a cell (identified by the 10X barcode), if there are 10K transcripts of one gene, each of the 10K transcripts will have a unique identifier. Thus, counting all unique identifiers per gene will gove us the absolute number of transcripts of that gene in that cell

The last / 4th part of the primer is a poly dT sequence which allows capture of a poly-A tailed RNA molecule in the cytoplasm or nucleus or any other compartment of the cell

The version 3 of 10X beads have additional sequence that allow capture of non poly A RNA molecules

The entire process happens on the Chromium machine, the GEMs are formed, after formation, the gel beads dissolve releasing the oligos into solution. Within the droplet, RNA is captured and converted to cDNA along with all their tags above. Then, the droplets are pooled together into bulk solution and a sequencing library is made. This library is then used to sequence on an Illumina machine. The data is then processed by 10X software and analysed.

10X recommends a sequencing depth of 20-50K reads per cell
