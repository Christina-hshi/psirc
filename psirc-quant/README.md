# psirc-quant
__psirc-quant__ is a program for quantifying abundances of both linear and circular full-length transcripts from RNA-seq data. It is part of the psirc project. __psirc-quant__ takes advantages of pseudo-alignment for efficient alignment of RNA-seq data to transcripts, and the likelihood maximization method for quantification in [kallisto](https://github.com/pachterlab/kallisto) which was designed to quantify only the linear transcripts.

## External libraries
- **zlib**
- **HDF5 C libraries**

## Installation
```
    git clone https://github.com/Christina-hshi/psirc-quant.git
    #go to project root directory
    cd psirc-quant
    mkdir release
    cd release
    cmake ..
    make psirc-quant
    #the psirc-quant program can be found at "src/psirc-quant"
```

## User manual
Run __psirc-quant__ without arguments for usage instructions. Following is the detailed manual to run __psirc-quant__ to quantify linear and circular transcripts in 2 steps.

We require the header lines of the circular transcripts in fasta format should end with "\tC" to let the program know that they are circular transcripts. And header lines of linear transcripts should not end with "\tC".

### Step 1: Build index
```
Usage: psirc-quant index [arguments] FASTA-files

Required argument:
-i, --index=STRING          Filename for the index to be constructed

Optional argument:
-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: 31)
    --make-unique           Replace repeated target names with unique names
```

### Step 2: Quantify
```
Usage: psirc-quant quant [arguments] FASTQ-files

Required arguments:
-i, --index=STRING            Filename for the index to be used for
                              quantification
-o, --output-dir=STRING       Directory to write output to

Optional arguments:
    --bias                    Perform sequence based bias correction
-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)
    --seed=INT                Seed for the bootstrap sampling (default: 42)
    --plaintext               Output plaintext instead of HDF5
    --fusion                  Search for fusions for Pizzly
    --single                  Quantify single-end reads
    --single-overhang         Include reads where unobserved rest of fragment is
                              predicted to lie outside a transcript
    --fr-stranded             Strand specific reads, first read forward
    --rf-stranded             Strand specific reads, first read reverse
-l, --fragment-length=DOUBLE  Estimated average fragment length
-s, --sd=DOUBLE               Estimated standard deviation of fragment length
-x, --min-fragment-length     Minimum length of a valid fragment
-X, --max-fragment-length     Maximum length of a valid fragment
                              (default: -l, -s values are estimated from paired
                               end data, but are required when using --single)
-t, --threads=INT             Number of threads to use (default: 1)
    --pseudobam               Save pseudoalignments to transcriptome to BAM file
    --genomebam               Project pseudoalignments to genome sorted BAM file
-g, --gtf                     GTF file for transcriptome information
                              (required for --genomebam)
-c, --chromosomes             Tab separated file with chrosome names and lengths
                              (optional for --genomebam, but recommended)
```

Authors
-------
- Christina Shi <hshi@cse.cuhk.edu.hk>
- Kevin Yip <kevinyip@cse.cuhk.edu.hk>
