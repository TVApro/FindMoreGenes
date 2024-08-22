import pandas as pd
import sys
from Bio import Entrez
import re

Entrez.email = 'A.N.Other@example.com'

#file = sys.argv[1]
file = 'Таблица.xlsx'
table = pd.ExcelFile(file)
df = table.parse('Sheet1')

for f in range(len(df.index)):
    with open('Списокгенов.txt', 'a') as output:
        print(df.loc[f]['Searched genes'], file=output)

with open('Списокгенов.txt', 'r') as inputs:
    for i in inputs.readlines():
        print(i)
        handle = Entrez.efetch(db="protein", id=i, rettype="fasta", retmode="text")
        sequence = handle.read()
        with open('Списокгенов.fasta', 'a') as output:
            print(sequence, file=output)

with open('Списокгенов.fasta', 'r') as inputs:
    for i in inputs.readlines():
        if i:
            with open('Списокгенов2.fasta', 'a') as output:                
                if '>' not in i:
                    print(i.replace('\n', ''), file=output)
                if '>' in i:
                    number = re.search(r'\w{2}_\d{9}.[0-9]|[A-Z0-9]{8}.[0-9]', i).group()
                    for f in range(len(df.index)):
                        if df.loc[f]['Searched genes'] == number:
                            name = df.loc[f]['Target organisms']
                    print('>'+number+' '+name)
                    print('>'+number+' '+name, file=output)
print('******************КОНЕЦ******************')
        
