********************************************************************************
MAST - Motif Alignment and Search Tool
********************************************************************************
	MAST version 3.5.1 (Release date: 2006/02/01 02:08:55)

	For further information on how to interpret these results or to get
	a copy of the MAST software please access http://meme.nbcr.net.
********************************************************************************


********************************************************************************
REFERENCE
********************************************************************************
	If you use this program in your research, please cite:

	Timothy L. Bailey and Michael Gribskov,
	"Combining evidence using p-values: application to sequence homology
	searches", Bioinformatics, 14(48-54), 1998.
********************************************************************************


********************************************************************************
DATABASE AND MOTIFS
********************************************************************************
	DATABASE -stdin (nucleotide)
	Database contains 4 sequences, 499297 residues

	Scores for positive and reverse complement strands are combined.

	MOTIFS /home/meme/TEST/tests/meme.adh.tcm (peptide)
	MOTIF WIDTH BEST POSSIBLE MATCH
	----- ----- -------------------
	  1    22   YCASKFAVRGFTRSMAMEYAPY
	  2    25   QGKVVLITGCSSGIGKATAKHFHKE

	PAIRWISE MOTIF CORRELATIONS:
	MOTIF     1
	----- -----
	   2   0.34
	No overly similar pairs (correlation > 0.60) found.

	Random model letter frequencies (from non-redundant database):
	A 0.070 C 0.024 D 0.040 E 0.052 F 0.040 G 0.074 H 0.029 I 0.041 K 0.052 
	L 0.096 M 0.017 N 0.032 P 0.065 Q 0.042 R 0.067 S 0.084 T 0.052 V 0.059 
	W 0.016 Y 0.022 
********************************************************************************


********************************************************************************
SECTION I: HIGH-SCORING SEQUENCES
********************************************************************************
	- Each of the following 4 sequences has E-value less than 10.
	- The E-value of a sequence is the expected number of sequences
	  in a random database of the same size that would match the motifs as
	  well as the sequence does and is equal to the combined p-value of the
	  sequence times the number of sequences in the database.
	- The combined p-value of a sequence measures the strength of the
	  match of the sequence to all the motifs and is calculated by
	    o finding the score of the single best match of each motif
	      to the sequence (best matches may overlap),
	    o calculating the sequence p-value of each score,
	    o forming the product of the p-values,
	    o taking the p-value of the product.
	- The sequence p-value of a score is defined as the
	  probability of a random sequence of the same length containing
	  some match with as good or better a score.
	- The score for the match of a position in a sequence to a motif
	  is computed by by summing the appropriate entry from each column of
	  the position-dependent scoring matrix that represents the motif.
	- Sequences shorter than one or more of the motifs are skipped.
	- The frame of the (best) motif match(es) is shown.
	  Frames 1, 2, and 3 are labeled a, b c, respectively.
	- The table is sorted by increasing E-value.
********************************************************************************

SEQUENCE NAME                      DESCRIPTION           FRAME   E-VALUE  LENGTH
-------------                      -----------           -----   -------- ------
gi|7296683|gb|AE003602.1|AE003602  Drosophila melanogaster... b    4.7e-18 297266
gi|7295475|gb|AE003567.1|AE003567  Drosophila melanogaster... c    4.3e-06 104808
gi|7289065|gb|AE002569.1|AE002569  Drosophila melanogaster... a       0.37  12850
gi|7289301|gb|AE002567.1|AE002567  Drosophila melanogaster... a          2  84373

********************************************************************************



********************************************************************************
SECTION II: MOTIF DIAGRAMS
********************************************************************************
	- The ordering and spacing of all non-overlapping motif occurrences
	  are shown for each high-scoring sequence listed in Section I.
	- A motif occurrence is defined as a position in the sequence whose
	  match to the motif has SEQUENCE p-value less than 0.0001.
	- The SEQUENCE p-value of a match is the probability of
	  some random subsequence in a set of n,
 where n is the sequence length minus the motif width plus 1,
	  scoring at least as well as the observed match.
	- For each sequence, all motif occurrences are shown unless there
	  are overlaps.  In that case, a motif occurrence is shown only if its
	  p-value is less than the product of the p-values of the other
	  (lower-numbered) motif occurrences that it overlaps.
	- The table also shows the E-value of each sequence.
	- Spacers and motif occurences are indicated by
	   o -d-    `d' residues separate the end of the preceding motif 
		    occurrence and the start of the following motif occurrence
	   o [snf]  occurrence of motif `n' with p-value less than 0.0001
		    in frame f.  Frames 1, 2, and 3 are labeled a, b c.
		    A minus sign indicates that the occurrence is on the
		    reverse complement strand.
********************************************************************************

SEQUENCE NAME                      E-VALUE   MOTIF DIAGRAM
-------------                      --------  -------------
gi|7296683|gb|AE003602.1|AE003602   4.7e-18  101659_[+2b]_1192_[+2c]_375_[+1c]_
                                             794_[-1b]_375_[-2b]_780_[-1b]_
                                             375_[-2b]_191218
gi|7295475|gb|AE003567.1|AE003567   4.3e-06  87179_[-1c]_17563
gi|7289065|gb|AE002569.1|AE002569      0.37  12850
gi|7289301|gb|AE002567.1|AE002567         2  84373

********************************************************************************



********************************************************************************
SECTION III: ANNOTATED SEQUENCES
********************************************************************************
	- The positions and p-values of the non-overlapping motif occurrences
	  are shown above the actual sequence for each of the high-scoring
	  sequences from Section I.
	- A motif occurrence is defined as a position in the sequence whose
	  match to the motif has SEQUENCE p-value less than 0.0001 as 
	  defined in Section II.
	- For each sequence, the first line specifies the name of the sequence.
	- The second (and possibly more) lines give a description of the 
	  sequence.
	- Following the description line(s) is a line giving the length, 
	  combined p-value, and E-value of the sequence as defined in Section I.
	- The next line reproduces the motif diagram from Section II.
	- The entire sequence is printed on the following lines.
	- Motif occurrences are indicated directly above their positions in the
	  sequence on lines showing
	   o the motif number/frame of the occurrence (a minus sign indicates that
	  the occurrence is on the reverse complement strand),
	   o the position p-value of the occurrence,
	   o the best possible match to the motif (or its reverse),
	   o columns whose match to the motif has a positive score (indicated 
	     by a plus sign), and
	   o the protein translation of the match (or its reverse).
********************************************************************************


gi|7296683|gb|AE003602.1|AE003602
  Drosophila melanogaster genomic scaffold 142000013386043 section 3 of 8, 
  complete sequence
  LENGTH = 297266  COMBINED P-VALUE = 1.17e-18  E-VALUE =  4.7e-18
  DIAGRAM: 101659_[+2b]_1192_[+2c]_375_[+1c]_794_[-1b]_375_[-2b]_780_[-1b]_
           375_[-2b]_191218


                                                  [+2b]
                                                  3.4e-13
                                                  Q..G..K..V..V..L..I..T..G..C..
                                                     +  +  +  +  +  +  +  +  +
                                                  A..G..K..V..V..L..I..T..G..A..
101617 CTTGAGATCCCCTTACATAATTCTGGATCGGATCATGAATTTCGCGGGCAAAGTGGTCCTTATTACGGGAGCA

       
       
       S..S..G..I..G..K..A..T..A..K..H..F..H..K..E..
       +  +  +  +  +     +  +  +  +     +  +  +  +
       S..S..G..I..G..A..A..T..A..I..K..F..A..K..Y..
101690 AGCTCCGGAATCGGAGCTGCAACCGCCATTAAGTTTGCCAAGTACGGCGCCTGTCTGGCTCTCAATGGACGCA

                                                                            [+2c
                                                                            9.8e
                                                                            Q..G
                                                                            +  +
                                                                            S..G
102858 AATGTAGTTTGGATTTTCATTTGAAATTTCGATTCGATGGGCAGCAAGGGCAAGGCTGGTCTGGACTTCTCCG

       ]
       -11
       ..K..V..V..L..I..T..G..C..S..S..G..I..G..K..A..T..A..K..H..F..H..K..E..
         +  +  +  +  +  +  +  +  +  +  +  +  +     +     +  +     +     +
       ..K..V..V..L..I..T..G..A..A..S..G..I..G..A..A..A..A..E..M..F..S..K..L..
102931 GCAAAGTGGTGCTTATCACGGGCGCAGCCTCCGGGATCGGGGCCGCCGCGGCGGAGATGTTCTCGAAGCTGGG

               [+1c]
               7.4e-08
               Y..C..A..S..K..F..A..V..R..G..F..T..R..S..M..A..M..E..Y..A..P..Y.
               +     +  +  +  +  +  +  +     +  +  +  +  +  +  +  +  +  +  +  +
               Y..N..M..S..K..A..A..V..D..Q..F..T..R..S..L..A..L..D..L..G..P..Q.
103369 TGGTGGCCTACAACATGTCCAAGGCGGCGGTGGACCAGTTTACCCGCTCCCTTGCCCTGGATCTGGGTCCCCA

       
       
       .
       
       .
103442 GGGTGTTCGGGTGAATGCGGTCAATCCGGGTGTGATTCGCACCAATCTGCAAAAGGCGGGGGGCATGGACGAG

                                                                        [-1b]
                                                                        7.4e-07
                                                                        ..Y..P..
                                                                          +  +
                                                                        ..K..P..
104172 AGTCCACCACGACGCTGCAGCTCGGTGATGATTACGCCGGGATTCACGGAGTTCACACGCACACCCTTGGGAG

       
       
       A..Y..E..M..A..M..S..R..T..F..G..R..V..A..F..K..S..A..C..Y
       +  +  +  +  +     +  +  +  +     +  +  +  +  +  +        +
       A..L..E..L..A..V..C..R..T..F..Q..D..V..A..A..K..S..V..N..Y
104245 CTAGCTCCAGAGCCACGCACCTGGTGAACTGATCCACGGCAGCCTTGGAAACATTGTATGCTAAGACTCCGGG

                                                                           [-2b]
                                                                           5.7e-
                                                                           ..E..
       
                                                                           ..L..
104610 GCCACTATCTGCTCCGCGGTCTCGTTGAGCTTATCCAAATTCCTGCCCACGATGGTGAGCAGGCCTCCCAGTT

       
       10
       K..H..F..H..K..A..T..A..K..G..I..G..S..S..C..G..T..I..L..V..V..K..G..Q
       +  +  +        +  +        +  +  +  +  +  +  +  +  +  +  +  +  +  +  +
       K..A..L..L..V..S..T..G..A..G..I..G..S..S..A..G..T..V..I..I..V..K..D..K
104683 TAGCCAAGAGCACCGAAGTACCCGCTCCAATTCCCGAACTGGCTCCGGTCACGATTATAACTTTATCCTTGAA

                                                      [-1b]
                                                      1.7e-05
                                                      ..Y..P..A..Y..E..M..A..M..
                                                        +  +  +  +  +  +  +
                                                      ..K..P..A..L..E..L..A..I..
105486 ATGTCAGTCACAATCACGCCGGGATTCACGGCGTTTACGCGGACTCCTTTGGGGGCCAGTTCCAGGGCTATGC

       
       
       S..R..T..F..G..R..V..A..F..K..S..A..C..Y
       +     +  +     +  +  +  +  +  +        +
       C..A..T..F..Q..D..V..A..A..K..S..V..N..Y
105559 AGGCCGTGAACTGGTCCACCGCTGCCTTAGACACATTGTAGGCCAGGACACCAGGAAAGGCACGCAGTCCACA

                                                         [-2b]
                                                         2.8e-11
                                                         ..E..K..H..F..H..K..A..
                                                              +  +  +  +     +
                                                         ..L..K..A..L..H..V..A..
105924 GTCTCCTTCAGCTTCTCCTCGTTGCGACCCACGATGACCAGGAGCCCTCCGAGTTTCGCCAAATGGACGGCGG

       
       
       T..A..K..G..I..G..S..S..C..G..T..I..L..V..V..K..G..Q
          +     +  +  +  +  +  +  +  +  +  +  +  +  +  +  +
       A..S..A..G..I..G..S..S..A..G..T..V..I..I..V..K..D..K
105997 CACTTGCTCCGATGCCGGAGCTGGCGCCGGTAACAATAATCACTTTGTCCTTGAAGGATGACATCGTGGGGTG


gi|7295475|gb|AE003567.1|AE003567
  Drosophila melanogaster genomic scaffold 142000013386050 section 54 of 54, 
  complete sequence
  LENGTH = 104808  COMBINED P-VALUE = 1.08e-06  E-VALUE =  4.3e-06
  DIAGRAM: 87179_[-1c]_17563


                        [-1c]
                        6.2e-08
                        ..Y..P..A..Y..E..M..A..M..S..R..T..F..G..R..V..A..F..K..
                          +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +     +
                        ..E..P..A..L..D..K..A..A..A..K..T..L..G..I..L..A..T..K..
87163  CAGTTGACGCGAATGCCCTCGGGCGCCAGATCCTTGGCGGCTGCCTTGGTCAAGCCAATCAGCGCGGTCTTGC

       
       
       S..A..C..Y
       +     +  +
       S..V..S..Y
87236  TGACGGAATAGGCTCCCAGTAGCTATGCGATTAAGATAACGGAGATAAGCATTGAACAACTGAACCGCAGATA


gi|7289065|gb|AE002569.1|AE002569
  Drosophila melanogaster genomic scaffold 142000013385354, complete sequence
  LENGTH = 12850  COMBINED P-VALUE = 9.28e-02  E-VALUE =     0.37
  DIAGRAM: 12850



gi|7289301|gb|AE002567.1|AE002567
  Drosophila melanogaster genomic scaffold 142000013385554, complete sequence
  LENGTH = 84373  COMBINED P-VALUE = 4.99e-01  E-VALUE =        2
  DIAGRAM: 84373


********************************************************************************


CPU: rocks-155.sdsc.edu
Time 1.352794 secs.

mast meme.adh.tcm -stdin -dna -seqp
