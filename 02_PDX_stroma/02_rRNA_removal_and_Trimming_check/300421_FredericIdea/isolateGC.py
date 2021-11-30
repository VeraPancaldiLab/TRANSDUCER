#!/usr/bin/env python
# coding: utf-8

# In[2]:


# imports
import os
import gzip
from Bio import SeqIO
from Bio.SeqUtils import GC


# In[7]:


file50 = 'PDAC221Tkf1_fwd_P50.fq'
file60 = 'PDAC221Tkf1_fwd_P60.fq'

# remove files in case they already exist so it dont get appended
try:
    os.remove(file50)
except:
    print(file50)
try:
    os.remove(file60)
except:
    print(file60)

# read fastq
with gzip.open('PDACtryGC.fq.gz', "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        readGC=GC(record.seq)
        
        if 46 < readGC < 51:
            with open(file50, 'a') as the_file:
                the_file.write(record.format("fastq"))
    
        elif 54 < readGC < 63:
            with open(file60, 'a') as the_file:
                the_file.write(record.format("fastq"))

