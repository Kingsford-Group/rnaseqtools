This repo `rnaseqtools` provides a set of tools to process transcripts (mainly in
`gtf` format).  To compile these tools, you first need to 
download the source code of the latest release
from [here](https://github.com/Kingsford-Group/rnaseqtools/releases/download/v1.0.2/rnaseqtools-1.0.2.tar.gz),
then use the following commands to compile:
```
./configure --prefix=/path/to/your/install/folder
make
make install
```

# gtfmerge
This tool is to compute the union of a collection of sets of transcripts.
Two transcripts are defined as identical if they are from the same chromosme,
on the same strand, and having the same intron chain coordinates.
The usage of `gtfmerge` is as follows:
```
gtfmerge union <input-gtf-list> <output-gtf-file> [-t <integer>] [-n]
```
1. The parameter `input-gtf-list` is mandatory, which provides a list of `gtf` files (each line specifies a file name).
Each `gtf` file gives a set of transcripts. 
2. The parameter `output-gtf-file` is mandatory, which contain the merged transcripts, also in `gtf` file format.
3. `-t <integer>` is optional. If it is provided, then the multiple-threading mode will be open, and the specified
number of threads will be used.
4. `-n` is optional. If it is used, then the number of appearance of each unioned transcript (i.e., how many input gtf files
contain this transcript) will be recorded and reported in the `cov` field of the output file. If this parameter is not used,
then the sum of the coverage of each unioned transcript (i.e., sum up of the coverage of all transcripts in the input
gtf files that are identical to this transcript) will be recorded and reported in the `cov` field of the output file.

# gtfcuff
This tool is to evaluate the accuracy of predicted transcripts. 
To use this tool, you first need to run `gffcompare`, which will
generate several files, and among them `gtfcuff` will usually
use `.tmap`. For example, to compute the AUC score (the parameter
used to draw the curve is the predicted coverage of all the transcripts),
you can use
```
gtfcuff auc <gffcompare.tmap> <number-of-exons-in-reference>
```
The last parameter is usually the number of multi-exon transcripts in the
reference annotation. You can find this number in the `.stats` file
produced by `gffcompare`. The AUC score will be printed to the standard output.
