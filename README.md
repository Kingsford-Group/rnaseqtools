This repo `rnaseqtools` provides a set of tools to process transcripts (mainly in
`gtf` format).  To compile these tools, you first need to clone
`scallop` repo and link the `lib` directory of `scallop` here. For example,
```
ln -sf ../scallop/lib .
```
After that you can run `./build.sh`, which shall generate the binary version 
of all tools.
