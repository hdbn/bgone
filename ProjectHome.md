# Be gone with the other arrays! #

The details are described in the following preprint: http://arxiv.org/abs/1310.1448




The runtime comparison for several algorithms, KKP1, 2, 3, LZISA6s, LZScan is put [here](https://docs.google.com/spreadsheet/ccc?key=0AlQjQSvkXsMmdG9land5ZmZnaWs1SXpqS2Q1QW9RS1E#gid=0).
(The competitor's implementations are [here](https://www.cs.helsinki.fi/group/pads/lz77.html) )

The comparison of peak memory consumptions of these algorithms is put [here](https://docs.google.com/spreadsheet/ccc?key=0AlQjQSvkXsMmdDlQR0JJVGFvS01FZ2JOUk5pNXREYWc#gid=0).

We only write the runtime of [divsufsort](https://code.google.com/p/libdivsufsort/) for the construction of Suffix Arrays, and all algorithms based on Suffix Arrays use divsufsort and includes the runtime of divsufsort in each runtimes.

Note, LZISA6s caused bus errors for some files, and we set lager value than 999 for the runtime of it in the case.
We set parameter of LZScan so that LZScan run to use memory of 6N bytes, 5N bytes, and 4N bytes, we call each runtime of them LZScan6, 5, 4 respectively.
Since it was hard to control the memory consumption accurately, real memory consumptions may be unexpected value.
Please see the comparison of peak memory consumptions for all algorithms.