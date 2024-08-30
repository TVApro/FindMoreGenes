import pandas as pd
import openpyxl
from io import StringIO
import re
import os
import sys
from bioservices import UniProt
from Bio import Entrez
import time

def zhest_start(i, idperc, e_val):
    print('Создание df из файла %s' % i)
    with open(i, 'r') as fileX:
        text = fileX.read()
        
    outfmt2 = text.find("Query=")
    if outfmt2 > 0:
        fmtsix = False
        print('It is not BLAST format 6. Aborting')
        exit()
    if outfmt2 < 0:
        fmtsix = True
        
    plus = lambda num: num + 1
    gene_list = text.split('\n')
    df = pd.DataFrame(columns=['Target genes', 'Target organisms', 'Searched genes', 'perc. ident (%)', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])

    if fmtsix == True:
         for i in gene_list:
            if i:
                str_list = i.split('\t')
                sch = 0
                for word in str_list:
                    if word:
                        sch = plus(sch)
                        if sch == 1:
                            qseqid = re.search(r'\|[A-aZ-z0-9]+\|', word).group().replace('|', '')
                        if sch == 2:
                            full_info = re.search(r'([A-aZ-z0-9_\-\.\\\/]+[\(A-aZ-z0-9_\-\.\\\)]*)', word).group()
                            number = re.search(r'REF[A-Z0-9]+|\w{2}_\d{9}.[0-9]|[A-Z0-9]{8}.[0-9]', word).group()
                            organism = full_info.replace(number, '').replace('_', ' ') 
                        if sch == 3:
                            pident = word
                        if sch == 4:
                            length = word
                        if sch == 5:
                            mismatch = word
                        if sch == 6:
                            gapopen = word
                        if sch == 7:
                            qstart = word
                        if sch == 8:
                            qend = word
                        if sch == 9:
                            sstart = word
                        if sch == 10:
                            send = word
                        if sch == 11:
                            evalue = word
                        if sch == 12:
                            bitscore = word
                            df.loc[len(df.index)] = [qseqid, organism, number, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore]                       
    def identity_persentage(YES, df, Title):
        drop_num = 0
        start_df_len = len(df.index)
        for f in range(len(df.index)):
            if float(df.loc[f][Title]) < float(YES):
                #print("%s %s %s dropped because of low identity persentage (%s)" % (f, df.loc[f]['Searched genes'], df.loc[f]['Target genes'], df.loc[f][Title]))
                df.drop(index=f, axis=0, inplace=True)
                drop_num += 1
        df = df.reset_index()
        df.drop("index", axis=1, inplace=True)
        print('\|/')
        print("%s of %s results was dropped" % (drop_num, start_df_len))
        return df

    def evalue(e_val, df, Title):
        drop_num = 0
        start_df_len = len(df.index)
        for f in range(len(df.index)):
            if df.loc[f][Title] != '0.0':
                grad = int(re.search('-[\d]+', df.loc[f][Title]).group().replace('-', ''))
                if grad <= int(e_val):
                    #print("%s %s %s dropped because of high evalue (%s)" % (f, df.loc[f]['Searched genes'], df.loc[f]['Target genes'], df.loc[f][Title]))
                    df.drop(index=f, axis=0, inplace=True)
                    drop_num += 1
        df = df.reset_index()
        df.drop("index", axis=1, inplace=True)
        print('\|/')
        print("%s of %s results was dropped" % (drop_num, start_df_len))
        return df

    if idperc.replace('.', '').isdigit() and float(idperc) <= 100 and float(idperc) > 0:
        print('<><><><><><><><><><><><><><><><><><><><><><><>')
        df = identity_persentage(idperc, df, 'perc. ident (%)')
    if e_val.isdigit() and int(e_val) >= 0:
        print('<><><><><><><><><><><><><><><><><><><><><><><>')
        df = evalue(e_val, df, 'evalue')
        
    return df

def zhest_finish(df):
    print('Анализ df')
    
    len1 = 0
    len2 = 1
    num = 0
    plus = lambda num: num + 1
    iteration = 0
    len0 = len(df.index)
    dupl_numb = 0
    while len1 != len2:
        len1 = len(df.index)
        num = plus(num)
        for h in range(len(df.index)):
            for g in range(len(df.index)):
                try:
                    #if df.loc[h]['Organism'] != "Methanothermobacter sp. K4": # unnessesary
                       #df.drop(index=h, axis=0, inplace=True)
                    if df.loc[h]['Searched genes'] == df.loc[g]['Searched genes'] and df.loc[h]['Target organisms'] == df.loc[g]['Target organisms']:
                        diff3 = float(df.loc[h]['perc. ident (%)']) - float(df.loc[g]['perc. ident (%)'])
                        diff2 = float(df.loc[h]['evalue']) - float(df.loc[g]['evalue'])
                        diff1 = float(df.loc[h]['bitscore']) - float(df.loc[g]['bitscore'])
                        diff_list = (diff1, diff2, diff3)
                        pass_num = 0
                        for diff in diff_list:
                            drop = False
                            while drop == False:
                                if diff > 0:
                                    #print("#%s | %s %s // %s ---> %s | %s " % (h, df.loc[g]['Target organisms'], df.loc[g]['Searched genes'], df.loc[g]['Target genes'], df.loc[h]['Target genes'], diff))
                                    df.drop(index=g, axis=0, inplace=True)
                                    drop = True
                                    break
                                if diff < 0:
                                    #print("#%s! | %s %s // %s ---> %s | %s " % (h, df.loc[h]['Target organisms'], df.loc[h]['Searched genes'], df.loc[h]['Target genes'], df.loc[g]['Target genes'], diff))
                                    df.drop(index=h, axis=0, inplace=True)
                                    drop = True
                                    break
                                if diff == 0:
                                    pass_num += 1
                                    if pass_num == 3:
                                        if df.loc[h]['Target genes'] != df.loc[g]['Target genes']:
                                            df.drop(index=g, axis=0, inplace=True)
                                        drop = True
                            break
                except:
                    pass
                
        #print('|||||||||||||||||||||||||||||||||||||||||||||')
        print('Iteration №%d is over' % num)
        len2 = len(df.index)
        print('There are %s results left of %s' % (len2, len1))
        #print('/////////////////////////////////////////////')
        df = df.reset_index()
        df.drop("index", axis=1, inplace=True)
    #print('Job is done')
    #print('There are %s results left of %s' % (len(df.index), len0))
    return df
   
def NCBI_what(x):
    ncbi_number = str(x)
    handle = Entrez.esummary(db="protein", id=ncbi_number)
    resultX = Entrez.read(handle)
    handle.close()
    return resultX[0]['Title']
    
def UNIPROT_what(x):
    u = UniProt(verbose=False)
    result = u.search(x, frmt='tsv', columns="gene_names, protein_name, organism_name")
    result_list = result.split('\n')
    result_result_list = re.split("\t", result_list[1])
    return result_result_list


# ЗДЕСЬ НАЧИНАЕТСЯ ТРЭШ
list_of_files = os.listdir(os.getcwd())
idperc = input('Set a minimum identity percentage? (No / 0.0000-100.0000): ')
e_val = input('Set a minimum e-value? (No / positive exponent of a fraction, like 70 in 3e-70): ')
YES = input("Get annotation data? (Yes, Да / Anything else): ")

file_name = sys.argv[1]+'_part_'
LISTOFNEEDED = list()
for file in list_of_files:
    if file_name in file:
        LISTOFNEEDED.append(file)
LISTOFNEEDED.sort()

start = True
for i in LISTOFNEEDED:
    if start == False:
        df_2 = zhest_finish(zhest_start(i, idperc, e_val))
        frames = [df, df_2]
        print('Объединяю датафреймы')
        df_3 = pd.concat(frames)
        df_3 = df_3.reset_index()
        df_3.drop("index", axis=1, inplace=True)
        df_3.drop_duplicates(inplace=True)
        if len(df_3.index) != 0:
            df = zhest_finish(df_3)
    if start == True:
        df = zhest_finish(zhest_start(i, idperc, e_val))
        start = False

yes_list = ('yes', 'y', 'Yes', 'Y', 'Д', 'Да','да', 'д', '', 'YES', 'ДА')
if YES not in yes_list:
    df.to_excel('Genes_of_interesting_without_annotation.xlsx')
    print("Done! Result was saved in Genes_of_interesting_without_annotation.xlsx")
    exit()

print()
print('<><><><><><><><><><><><><><><><><><><><><><><>')
print("Checking the names of genes in Uniprot and NCBI")
    
df.insert(loc=0, column='Host organisms', value="")
df["UNIPROT_query_name"] = ""
df["UNIPROT_prot_name"] = ""
df["NCBI_prot_name"] = ""
Entrez.email = 'A.N.Other@example.com'

for i in range(len(df.index)):
    try:
            UNIPROT_result = UNIPROT_what('"'+df.loc[i]['Target genes']+'"')
            UNIPROT_gene_name = UNIPROT_result[0]
            UNIPROT_protein_name = UNIPROT_result[1]
            UNIPROT_host_name = UNIPROT_result[2]
    except:
            print('Sorry, I was unable to access Uniprot :( ')
            UNIPROT_gene_name = "ND"
            UNIPROT_protein_name = "ND"
    try_number = 0
    try:
        try_number += 1
        NCBI_name = NCBI_what(df.loc[i]['Searched genes'].replace(' ', ''))
        print("%s | %s %s --> %s / %s" % (i, df.loc[i]['Target organisms'], df.loc[i]['Target genes'], UNIPROT_gene_name, NCBI_name))
    except:
        print('Oops, I was unable to get the NCBI number, trying №%s' % try_number)
        NCBI_name = '!'
        while NCBI_name == '!':
            try:
                    NCBI_name = NCBI_what(df.loc[i]['Searched genes'].replace(' ', ''))
                    print("%s | %s %s --> %s / %s" % (i, df.loc[i]['Target organisms'], df.loc[i]['Target genes'], UNIPROT_gene_name, NCBI_name))
            except:
                    try_number += 1
                    print('Oops, I was unable to get the NCBI number, trying №%s' % try_number)
                    if try_number == 5:
                            NCBI_name = "Sorry, I couldn't :( Don't angry, please!"
                            NCBI_name = "ND"
                            print("%s | %s %s --> %s / %s" % (i, df.loc[i]['Target organisms'], df.loc[i]['Target genes'], UNIPROT_gene_name, NCBI_name))

    df.at[i, 'Host organisms'] = UNIPROT_host_name
    df.at[i, 'UNIPROT_query_name'] = UNIPROT_gene_name
    df.at[i, "UNIPROT_prot_name"] = UNIPROT_protein_name
    df.at[i, "NCBI_prot_name"] = NCBI_name

    
df.to_excel('Genes_of_interesting.xlsx')
print("Result was saved in 'Genes_of_interesting.xlsx'")
print()
