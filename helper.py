import pandas as pd
import sys
from Bio import Entrez
from bioservices import UniProt
import re

#file = sys.argv[1]
file = 'PFL-like.txt'
u = UniProt(verbose=False)
Entrez.email = 'A.N.Other@example.com'

def def_1(file):
    new_set = set()
    with open(file, 'r') as inputs_1:
        for i in inputs_1.readlines():
            if i:
                new_set.add(i.replace('\n', ''))
    for j in new_set:
        print(j)
        handle = Entrez.efetch(db="protein", id=j, retmode='xml')
        LOCUS = Entrez.read(handle)[0]["GBSeq_locus"]
        handle.close
        with open('INPUT_LOCUS.txt', 'a') as inputs_2:
            print(LOCUS, file=inputs_2)
def def_2():
    with open('INPUT_LOCUS.txt', 'r') as inputs_3:
        for LOCUS in inputs_3.readlines():
            if LOCUS:
                fasta = u.search(LOCUS, frmt='fasta')
                print(fasta)
                with open('INPUT.fasta', 'a') as inputs_4:
                    print(fasta, file=inputs_4)
        
def_1(file)
def_2()
print('******************КОНЕЦ******************')
        
