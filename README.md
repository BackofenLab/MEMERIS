[![Build Status](https://travis-ci.org/BackofenLab/MEMERIS.svg?branch=master)](https://travis-ci.org/BackofenLab/MEMERIS)
# MEMERIS
MEMERIS is an extension of the MEME (Bailey and Elkan (1994)) motif finder. 
MEMERIS integrates information about RNA secondary structures into the motif search to guide the search towards single-stranded regions. 

**Abstract**: RNA binding proteins recognize RNA targets in a sequence specific manner. Apart from the sequence, the secondary structure context of the binding site also affects the binding affinity. Binding sites are often located in single-stranded RNA regions and it was shown that the sequestration of a binding motif in a double-strand abolishes protein binding. Thus, it is desirable to include knowledge about RNA secondary structures when searching for the binding motif of a protein. We present the approach MEMERIS for searching sequence motifs in a set of RNA sequences and simultaneously integrating information about secondary structures. To abstract from specific structural elements, we precompute position-specific values measuring the single-strandedness of all substrings of an RNA sequence. These values are used as prior knowledge about the motif starts to guide the motif search. Extensive tests with artificial and biological data demonstrate that MEMERIS is able to identify motifs in single-stranded regions even if a stronger motif located in double-strand parts exists. The discovered motif occurrences in biological datasets mostly coincide with known protein-binding sites. This algorithm can be used for finding the binding motif of single-stranded RNA-binding proteins in SELEX or other biological sequence data.

## Contribution

Feel free to contribute to this project by wirting [Issues](https://github.com/BackofenLab/MEMERIS/issues) with feature requests or bug reports.

## Cite
If you use MEMERIS, please cite our [article](http://nar.oxfordjournals.org/content/34/17/e117):
```
doi: 10.1093/nar/gkl544
```
