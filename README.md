===============
About Bamfilter
===============

Bamfilter is a fast filtering tool for BAM files.
Generally, it works as a part of ReALLEN package.

============
Requirements
============

BamTools API is required.  
https://github.com/pezmaster31/bamtools/

A quick way is to add the included library path to the environment variable 'LD_LIBRARY_PATH'.

    export LD_LIBRARY_PATH=<install dir>/include/bamtools/lib:$LD_LIBRARY_PATH

============
Installation
============

    git clone https://github.com/cannaryo/bamfilter.git
    cd bamfilter
    make
    make install

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

========
Citation
========

* Bamfiler is developped for faster filtering in ReALLEN workflow.  
ReALLEN: structural variation discovery in cancer genome by sensitive analysis of single-end reads.  
 Ryo Kanno, Daisuke Tanaka, Hideaki Nanamiya, Takao Isogai
