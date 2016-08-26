#!/usr/local/bin/python
import sys, linecache
count = linecache.getline(sys.argv[1],int (sys.argv[2]))
outf = open('out.dat', 'w')          # open for output (creates)
outf.write(count)        # write a line of text
outf.close()
