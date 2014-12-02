PGA
===

So far this repository contains code to produce a hilbert curve based representation of genomic variants identified by GATK.
First the GATK output needs to be parsed (parseVCF.sh) to put the variant information in a format needed for the visualization 
code (GenomeViz.R). So far GENE_NAME cannot be included but this is a TODO.
The visualizing code (GenomeViz.R) assumes that about 4 Million high quality (PASS) variants are reported in the VCF file. If 
significantly more or less variants are reported the hilbert curve level (standard is 11) needs to be increased (to 12) or 
decreased (to 10) in the code (GenomeViz.R).

The output is a series of .png files plotting the variants in hilbert representation in different zoom stages (corresponding to 
hilbert levels 0-11). Each zoom stage is represented with and without focus rectangle.
Additionally the legend explaining the color code and stating the number of variants in each category is plotted in 
Sample_legend.png.

The .png files can be loaded in an animated gif maker (e.g. http://sourceforge.net/projects/gifapp/) and put together as 
an animated gif. 1 second per frame has been found to be a good speed.    
