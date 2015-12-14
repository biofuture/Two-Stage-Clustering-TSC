# Two-Stage-Clustering-TSC
This C++ program is used to cluster 16S rRNA sequences into clusters (OTUs) with a two stage strategy, by sperating the sequences into
abundant and rare by a cut off number, the distance between pair of sequences is calculated by a kmer based reduction method reduce the
computation time. The pair-wise global distance is needleman wunsch algorithm using the package "seqan". So If users try to compile the
program from source code, Users need to download the "seqan" package (http://packages.seqan.de/). The lower abundant sequences are 
processed by a greey algorithm to reduce the computational time.

If you use this program, please cite: 
    Jiang, X.-T., H. Zhang, H.-F. Sheng, Y. Wang, Y. He, F. Zou and H.-W. Zhou (2012). "Two-stage clustering (TSC): a pipeline for selecting operational taxonomic units for the high-throughput sequencing of PCR amplicons." PLoS One 7(1): e30230-e30230.

        Two Stage Clustering Version 1.2 CMS by GIT
        Last Modified 22/11/2011
        by Jiang Xiao-Tao
        biofuture.jiang@gmail.com/hust.cn@163.com
        Usage: ./TSC -i <.fasta> -o <out.prefix> -c [abundancecutoff]
        -x [terminal gap] -n [numberofthread] -k [ksize] -f [kmerfilter]
        -s [NotMatchCut] -d [distance] -m [al/cl/sl]

        -i  <str>   input fasta file
        -o  <str>   output file prefix
        -c [number] the cutoff for high and low abundant default 3
        -g [gapopen] gapopen penality default -10
        -e [gapextension] default -1
        -a [match] substiution matrix match default 1
        -b [mismatch] substitution matrix mismatch default -1
        -k [number] the word size for kmer default 6
        -f [double] the kmer filter default 0.5
        -d [double] the distance to do OTU picking default 0.03
        -n [number] number of thread used default 4
        -x [0/1] for 454 data if 1 then caculate end gap         if 0 remove end gap default 1
        -s [number] the kind number of the kmerdist not fit the distance require default 10
        -m [method] the clust method for high abundant sequences
        -r [data type] the data type to process deafult illumina data if using 454 data set -r 454

Have fun !
  
