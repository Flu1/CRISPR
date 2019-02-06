# CRISPR
Code used to identify successful CRISPR clones for Staller et al. 2019.

Running this code in R should allow you to identify potential clones which have successfully had genes knocked out by CRISPR.
The full sequencing files are available from the ENA under project PRJEB31093 and the barcodes used can be found in the attached spreadsheet.

To run the code, you will need to run the program barcodecco in R. I would recommend testing with pool 21.  Download the sequences from the ENA (https://www.ebi.ac.uk/ena) project reference PRJEB31093 and map them to a reference of ANP32A and B (e.g. use map to reference in Geneious). The reference file is called "ANP32A and B ref.fasta". Save the file as a sam file called pool21.sam and then run the following command barcodecco(pool21.sam,"CGATTG"). 
This should download all the dependencies and give the two alleles for this barcode that were above the cutoff.
