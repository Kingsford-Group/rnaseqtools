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
