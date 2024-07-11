NCBI-Bioprojects
PRJDB4868
PRJNA290784
PRJNA259658

to download the data:
go to:
https://www.ncbi.nlm.nih.gov/datasets/genome/?bioproject=PRJDB4868
https://www.ncbi.nlm.nih.gov/datasets/genome/?bioproject=PRJNA290784
https://www.ncbi.nlm.nih.gov/datasets/genome/?bioproject=PRJNA259658

then, click the box next to "Assembly" in the first column to select genomes.
Then download packages for RefSeq genomes, and only select the box
to download the GBFF sequence and annotation files.

Download each package into the data/PRJDB4868, data/PRJNA290784, and data/PRJNA259658
directories, respectively. Then, put the directories containing the genomic.gbff files
(each of these is named by RefSeq "GCF_XXXXX" IDs)
directly into the data/PRJDB4868 directory for example, removing the extra layers of
subdirectories.

If this is confusing, check the python code to see how the directory structures are parsed.