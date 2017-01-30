===============
About Bamfilter
===============

Bamfilter is a fast filtering tool for BAM files.
Generally, it works as a part of ReALLEN package.

============
Installation
============

    git clone https://github.com/cannaryo/bamfilter.git
    cd bamfilter
    make
    make install

============
Requirements
============

BamTools API is required.  
A quick way is to add the library path under Bamfilter directory to the environment variable *LD_LIBRARY_PATH*.

    export LD_LIBRARY_PATH=<bamfilter-root-dir>/include/bamtools/lib:$LD_LIBRARY_PATH

Or, compile the library from original sources.  
https://github.com/pezmaster31/bamtools/

==========================
Run Bamfilter with ReALLEN
==========================    

    reallen-fast.sh <your-bam-file-name.bam>

You must have installed ReALLEN package to use this command.  
See, 
https://github.com/cannaryo/reallen/

==========================
Run Bamfilter individually
==========================

See help messages for detailed information.

    bamfilter --help

Example of the simple filtering:

    bamfilter <your-bam-file-name.bam> -u -p -b 8 -s 8 -o <output-file.sam>

Example of the operation with soft-clip split and fixed-length split:

    bamfilter <your-bam-file-name.bam> -u -p -b 8 -s 8 -k -o <output-file.um.sam> \
              --softclip <output-file.sc.fq> \
              --fixed <output-file.um.fq> -f 45 -l 100

*--softclip* option enables soft-clip split.  
*--fixed* opttion enables fixed-length split.

