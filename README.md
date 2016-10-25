# **PUNASfilter**

We present a parallel ungapped-alignment-featured seed verification (PUNAS) algorithm,
a fast filter for effectively removing the majority of false positive seeds, thus significantly 
accelerating the short-read alignment process. PUNAS is based on bit-parallelism and takes
advantage of SIMD vector units of modern microprocessors.  PUNAS algorithm expands the basic Hamming 
distance computation, which only detects substitutions, into a full-fledged edit-distance filter, which 
counts not only substitutions but insertions and deletions as well. 

We have implemented our PUNAS algorithm as a standalone application (PUNASfilter) that can be 
incorporated into existing seed-and-extend based SRA tools. PUNASfilter is a filter that detects and 
filters some of the string pairs that have edit-distances that are greater than error threshold T, but it 
does not validate the string pairs that pass the filter regarding if they have edit-distances smaller than 
T. In a DNA mapper or protein mapper, PUNASfilter should then be extended by a more expensive DP-based  
alignment method. The PUNASfilter application is written in OpenMPaugmented C++. Our implementation 
employs a vectorize-and-scale approach and can be execyted on multi-core CPUs and many-core Knights 
Landing (KNL)-based Xeon Phi processors as well.

Currently, we presents three implementations: they are PUNAS-SSE, PUNAS-AVX2, PUNAS-AVX512. They are 
implemented using Intel SSE, AVX2 and AVX512 instructuon set respectively. Besides, the query string 
length is not limited by vector registers size. Although currently PUNAS only supports mapping DNA, it can 
be easily expanded to support protein matching. Please directly contact the author if you need such 
functionality. 

# **Getting started**

To build PUNAS, after you choose one in three implementation, based on your hardware platform, then simply 

do:
```bash
$make
```
# **Running a test**

To run a test using PUNAS, simply do:

```bash
$ ./PUNAS <reference file> <read file> <error threashold>
```
After the execution is finished, a summary of the run is printed out, including the total number of string 
pairs processed and the total number of string pairs pass PUNAS.



