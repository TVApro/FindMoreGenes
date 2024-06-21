import pandas as pd
import openpyxl
from io import StringIO
import re

myfile = "PSIBLAST-OUTPUT.txt"

with open(myfile, 'r') as fileX:
    text = fileX.read()

outfmt2 = text.find("Query=")
if outfmt2 > 0:
    print('OUTFMT #2')
    fmtsix = False
if outfmt2 < 0:
    print('OUTFMT #6')
    fmtsix = True

stop = lambda num: num + 1
num = 0
len1 = 0
len2 = 1
iteration = 0
    
if fmtsix == False:
    gene_list = text.split("Query=")
    df = pd.DataFrame(columns=['Query', 'Query_number', 'Query_org', 'Organism', 'Gene_number', 'Bitscore', 'Evalue'])
if fmtsix == True:
    gene_list = text.split('\n')
    df = pd.DataFrame(columns=['qseqid', 'Organism', 'Gene_number', 'Perc. ident (%)', 'lenth', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])

if fmtsix == False:    
    for i in gene_list:
        if i:
            try:
                zag = re.search(r'OS=\w+ \w+', i.replace('\n', ' ')).group().replace('OS=', '')
                gen = re.search(r'GN=\w+', i.replace('\n', ' ')).group().replace('GN=', '')
                gen_num = re.search(r'tr\|[A-Za-z0-9]+\|[A-Za-z0-9]+_[A-Za-z0-9]+', i.replace('\n', ' ')).group()
                i_list = i.split('\n')
                for j in i_list:
                    if j:
                        try:
                            number = re.search(r'WP_\d{9}.\d{1}', j).group()
                            exist=True
                        except:
                            exist=False
                            pass
                        if exist == True:
                            organism = re.search(r'(\w+_)+(sp._[A-aZ-z0-9]*_WP_)|(\w+_)(\w*_WP_)', j).group().replace('_WP_', '').replace('_', ' ')
                            bitscore = re.search(r'\s+\d+(.\d)?\s+', j).group().replace(' ', '')
                            evalue = re.search(r'\d+(.\d)?(\s{4})(\s{1})?\d+(.\d)?(e-\d+)?', j).group()[4:].replace(' ', '')
                            df.loc[len(df.index)] = [gen, gen_num, zag, organism, number, bitscore, evalue]
            except:
                pass
if fmtsix == True:
     for i in gene_list:
        if i:
            str_list = i.split('\t')
            sch = 0
            for word in str_list:
                if word:
                    sch = stop(sch)
                    if sch == 1:
                        qseqid = word
                    if sch == 2:
                        organism = re.search(r'([A-aZ-z0-9_\.\\\/]+_WP_)', word).group().replace('_WP_', '').replace('_', ' ')
                        number = re.search(r'WP_\d{9}.\d{1}', word).group()
                    if sch == 3:
                        pident = word
                    if sch == 4:
                        lenth = word
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
                        df.loc[len(df.index)] = [qseqid, organism, number, pident, lenth, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore]                       
len0 = len(df.index)
dupl_numb = 0

while len1 != len2:
    len1 = len(df.index)
    num = stop(num)
    for h in range(len(df.index)):
        for g in range(len(df.index)):
            try:
                if fmtsix == False:
                    if df.loc[h]['Gene_number'] == df.loc[g]['Gene_number']:
                        if df.loc[h]['Query'] != df.loc[g]['Query'] or df.loc[h]['Query_org'] != df.loc[g]['Query_org'] or df.loc[h]['Query_number'] != df.loc[g]['Query_number']:
                            diff = float(df.loc[h]['Bitscore']) - float(df.loc[g]['Bitscore'])
                            diff2 = float(df.loc[h]['Evalue']) - float(df.loc[g]['Evalue'])
                            if diff > 0:
                                print("The result %s %s %s was rejected" % (df.loc[g]['Query'], df.loc[g]['Gene_number'], df.loc[g]['Bitscore']))
                                print("    in favor of the %s %s %s " % (df.loc[h]['Query'], df.loc[h]['Gene_number'], df.loc[h]['Bitscore']))   
                                df.drop(index=g, axis=0, inplace=True)
                            if diff < 0:
                                print("The result %s %s %s was rejected" % (df.loc[h]['Query'], df.loc[h]['Gene_number'], df.loc[h]['Bitscore']))
                                print("    in favor of the %s %s %s " % (df.loc[g]['Query'], df.loc[g]['Gene_number'], df.loc[g]['Bitscore']))  
                                df.drop(index=h, axis=0, inplace=True)
                            if diff == 0:
                                print('bitscores are equal')
                                diff2 = float(df.loc[g]['Evalue']) - float(df.loc[h]['Evalue'])
                                if diff2 > 0:
                                    print("The result %s %s %s was rejected" % (df.loc[g]['Query'], df.loc[g]['Gene_number'], df.loc[g]['Evalue']))
                                    print("    in favor of the %s %s %s because of evalue" % (df.loc[h]['Query'], df.loc[h]['Gene_number'], df.loc[h]['Evalue']))   
                                    df.drop(index=g, axis=0, inplace=True)
                                if diff2 < 0:
                                    print("The result %s %s %s was rejected" % (df.loc[h]['Query'], df.loc[h]['Gene_number'], df.loc[h]['Evalue']))
                                    print("    in favor of the %s %s %s because of evalue" % (df.loc[g]['Query'], df.loc[g]['Gene_number'], df.loc[g]['Evalue']))  
                                    df.drop(index=h, axis=0, inplace=True)
                                if diff2 == 0:
                                    print("evalues are equal")
            except:
                pass
            try:
                if fmtsix == True:
                    #if df.loc[h]['Organism'] != "Methanothermobacter sp. K4": # unnessesary
                       #df.drop(index=h, axis=0, inplace=True)
                    if df.loc[h]['Gene_number'] == df.loc[g]['Gene_number'] and df.loc[h]['Organism'] == df.loc[g]['Organism']:
                        diff3 = float(df.loc[h]['Perc. ident (%)']) - float(df.loc[g]['Perc. ident (%)'])
                        diff2 = float(df.loc[h]['evalue']) - float(df.loc[g]['evalue'])
                        diff1 = float(df.loc[h]['bitscore']) - float(df.loc[g]['bitscore'])
                        diff_list = (diff1, diff2, diff3)
                        pass_num = 0
                        for diff in diff_list:
                            drop = False
                            while drop == False:
                                if diff > 0:
                                    print("#%s | %s %s // %s ---> %s | %s " % (h, df.loc[g]['Organism'], df.loc[g]['Gene_number'], df.loc[g]['qseqid'], df.loc[h]['qseqid'], diff))
                                    df.drop(index=g, axis=0, inplace=True)
                                    drop = True
                                    break
                                if diff < 0:
                                    print("#%s! | %s %s // %s ---> %s | %s " % (h, df.loc[h]['Organism'], df.loc[h]['Gene_number'], df.loc[h]['qseqid'], df.loc[g]['qseqid'], diff))
                                    df.drop(index=h, axis=0, inplace=True)
                                    drop = True
                                    break
                                if diff == 0:
                                    pass_num += 1
                                    if pass_num == 3:
                                        drop = True
                            break
            except:
                pass
            
    print('|||||||||||||||||||||||||||||||||||||||||||||')
    print('Iteration №%d is over' % num)
    len2 = len(df.index)
    print('There are %s results left of %s' % (len2, len1))
    print('/////////////////////////////////////////////')
    df = df.reset_index()
    df.drop("index", axis=1, inplace=True)
print('Job is done')
print('There are %s results left of %s' % (len(df.index), len0))

def identity_persentage(YES, df, Title):
    drop_num = 0
    for f in range(len(df.index)):
        if float(df.loc[f][Title]) < float(YES):
            print("%s %s %s dropped because of low identity persentage (%s)" % (f, df.loc[f]['Gene_number'], df.loc[f]['qseqid'], df.loc[f][Title]))
            df.drop(index=f, axis=0, inplace=True)
            drop_num += 1
    df = df.reset_index()
    df.drop("index", axis=1, inplace=True)
    print('\|/')
    print("%s of %s results was dropped" % (drop_num, len(df.index)))
    return df

def evalue(e_val, df, Title):
    drop_num = 0
    for f in range(len(df.index)):
        if df.loc[f][Title] != '0.0':
            grad = int(re.search('-[\d]+', df.loc[f][Title]).group().replace('-', ''))
            if grad <= int(e_val):
                print("%s %s %s dropped because of high evalue (%s)" % (f, df.loc[f]['Gene_number'], df.loc[f]['qseqid'], df.loc[f][Title]))
                df.drop(index=f, axis=0, inplace=True)
                drop_num += 1
    df = df.reset_index()
    df.drop("index", axis=1, inplace=True)
    print('\|/')
    print("%s of %s results was dropped" % (drop_num, len(df.index)))
    return df

if fmtsix == True:
    print()
    idperc = input('Set a minimum identity percentage? (No / 0.0000-100.0000): ')
    if idperc.replace('.', '').isdigit() and float(idperc) <= 100 and float(idperc) > 0:
        print('<><><><><><><><><><><><><><><><><><><><><><><>')
        df = identity_persentage(idperc, df, 'Perc. ident (%)')
        print('Inappropriate number. Pass')

    print()
    e_val = input('Set a minimum e-value? (No / positive exponent of a fraction, like 70 in 3e-70): ')
    if e_val.isdigit() and int(e_val) >= 0:
        print('<><><><><><><><><><><><><><><><><><><><><><><>')
        df = evalue(e_val, df, 'evalue')
    else:
        print('Inappropriate number. Pass')

    print()
    YES = input("Get annotation data? (Yes)")
    yes_list = ('yes', 'y', 'Yes', 'Y', 'Д', 'Да','да', 'д', '', 'YES', 'ДА')
    if YES not in yes_list:
        df.to_excel('Genes_of_interesting_6f_without_annotation.xlsx')
        exit()
    print()
    print('<><><><><><><><><><><><><><><><><><><><><><><>')
    print("Checking the names of genes in Uniprot and NCBI")
    from bioservices import UniProt
    from Bio import Entrez
    import time
    
    def NCBI_what(x):
    	ncbi_number = str(x)
    	handle = Entrez.esummary(db="protein", id=ncbi_number)
    	resultX = Entrez.read(handle)
    	handle.close()
    	return resultX[0]['Title']
        
    def UNIPROT_what(x):
    	u = UniProt(verbose=False)
    	result = u.search(x, frmt='tsv', columns="gene_names, protein_name")
    	result_list = result.split('\n')
    	result_result_list = re.split("\t", result_list[1])
    	return result_result_list
        
    df["UNIPROT_query_name"] = ""
    df["UNIPROT_prot_name"] = ""
    df["NCBI_prot_name"] = ""
    Entrez.email = 'A.N.Other@example.com'
    
    for i in range(len(df.index)):
        try:
        	UNIPROT_result = UNIPROT_what('"'+re.search(r'\|[A-aZ-z0-9]+\|', df.loc[i]['qseqid']).group().replace('|', '')+'"')
        	UNIPROT_gene_name = UNIPROT_result[0]
        	UNIPROT_protein_name = UNIPROT_result[1]
        except:
       		print('Sorry, I was unable to access Uniprot :( ')
       		UNIPROT_gene_name = "ND"
       		UNIPROT_protein_name = "ND"
        try_number = 0
        try:
            try_number += 1
            NCBI_name = NCBI_what(df.loc[i]['Gene_number'].replace(' ', ''))
            print("%s | %s %s --> %s / %s" % (i, df.loc[i]['Organism'], df.loc[i]['qseqid'], UNIPROT_gene_name, NCBI_name))
        except:
            print('Oops, I was unable to get the NCBI number, trying №%s' % try_number)
            NCBI_name = '!'
            while NCBI_name == '!':
            	try:
            		NCBI_name = NCBI_what(df.loc[i]['Gene_number'].replace(' ', ''))
            		print("%s | %s %s --> %s / %s" % (i, df.loc[i]['Organism'], df.loc[i]['qseqid'], UNIPROT_gene_name, NCBI_name))
            	except:
            		try_number += 1
            		print('Oops, I was unable to get the NCBI number, trying №%s' % try_number)
            		if try_number == 5:
            			NCBI_name = "Sorry, I couldn't :( Don't angry, please!"
            			NCBI_name = "ND"
            			print("%s | %s %s --> %s / %s" % (i, df.loc[i]['Organism'], df.loc[i]['qseqid'], UNIPROT_gene_name, NCBI_name))
 
        df.at[i, 'UNIPROT_query_name'] = UNIPROT_gene_name
        df.at[i, "UNIPROT_prot_name"] = UNIPROT_protein_name
        df.at[i, "NCBI_prot_name"] = NCBI_name
        
    df.to_excel('Genes_of_interesting_6f.xlsx')
    print("Done! Result was saved in Genes_of_interesting_6f.xlsx")


if fmtsix == False:
    df.to_excel('Genes_of_interesting_2f.xlsx')
    print("Done! Result was saved in Genes_of_interesting_2f.xlsx")
