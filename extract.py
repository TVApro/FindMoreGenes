import pandas as pd
import sys
from Bio import Entrez
from bioservices import UniProt
import re

Entrez.email = 'A.N.Other@example.com'

#file = sys.argv[1]
file = 'Genes_of_interesting.xlsx'
table = pd.ExcelFile(file)
df = table.parse('Sheet1')

def def_1(df):
    for f in range(len(df.index)):
        with open('List_of_genes_1.txt', 'a') as output:
            print(df.loc[f]['Searched genes'], file=output)
    print("---> Файл с номерами белков 'List_of_genes_1.txt' записан")
    
    
def def_2():
    print('Получение последовательностей:')
    with open('List_of_genes_1.txt', 'r') as inputs:
        for i in inputs.readlines():
            print(i.replace('\n', ''))
            if 'REF' not in i:
                handle = Entrez.efetch(db="protein", id=i, rettype="fasta", retmode="text")
                sequence = handle.read()
                with open('List_of_genes_1.fasta', 'a') as output:
                    print(sequence, file=output)
            else:
                u = UniProt(verbose=False)
                sequence = u.search(i.replace('REF', ''), frmt='fasta')
                with open('List_of_genes_1.fasta', 'a') as output:
                    print(sequence, file=output)
    print()
    print("---> Файл с последовательностями белков 'List_of_genes_1.fasta' записан")
    
def def_3():
    print("Упрощение имён:")
    ref_list = {}
    with open('List_of_genes_1.fasta', 'r') as inputs:
        for i in inputs.readlines():
            if i:
                with open('List_of_genes_2.fasta', 'a') as output:                
                    if '>' not in i:
                        print(i.replace('\n', ''), file=output)
                    if '>' in i:
                        try:
                            number = re.search(r'\w{2}_\d{9}.[0-9]|[A-Z0-9]{8}.[0-9]', i).group()
                        except:
                            number = 'REF_'+re.search(r'\|[A-Z0-9]+\|', i).group().replace('tr', '').replace('|', '')
                        for f in range(len(df.index)):
                            if df.loc[f]['Searched genes'] == number:
                                name = df.loc[f]['Target organisms']
                                length = df.loc[f]['length']
                                ref_name = df.loc[f]['Host organisms']
                                ref_gene = df.loc[f]['Target genes']
                                ref_list[ref_name] = ref_gene
                        name_2 = '>'+number+' '+name
                        print(name_2.replace(' ', '_'))
                        print(name_2.replace(' ', '_'), file=output)
    print()
    print("---> Записан файл FASTA 'List_of_genes_2.fasta' с модифицированными названиями")
    
    
    return ref_list

def def_4(ref_list):
    print("Поиск референсов:")
    with open('List_of_genes_2.fasta', 'a') as output: 
        for j in ref_list:
            u = UniProt(verbose=False)
            result = u.search(ref_list[j].replace('REF_', ''), frmt='fasta')
            result_fasta = ''.join(result.split('\n')[1:])
            name = '>REF'+ref_list[j]+' '+j
            print(name.replace(' ', '_'))
            print(name.replace(' ', '_')+'_'+str(len(result_fasta))+'_aa', file=output)
            print(result_fasta, file=output)
    print("---> В файл FASTA 'List_of_genes_2.fasta' добавлены референсные последовательности из UniProt")
    print()

print('Извлекаю последовательности белков из результирующего файла')
print()
def_1(df)
def_2()
ref_list = def_3()
def_4(ref_list)
print('******************КОНЕЦ******************')
        
