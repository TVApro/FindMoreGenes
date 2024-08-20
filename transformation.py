import os
import re
import shutil

def fasta_for_analisys(text, G_name):
    fasta_list = text.split(">")
    FASTA ={}
    for i in fasta_list:
        if i:
            end = i.find("\n")
            prot_number = re.search(r'\w{2}_\d{9}.[0-9]|[A-Z0-9]{8}.[0-9]', i[:end]).group() # to find number of protein
            fasta_name = '>' + G_name + '_'+prot_number
            fasta_sequence = re.sub('\n', '', i[end:])
            FASTA[fasta_name]=fasta_sequence
    return FASTA

def print_fasta(FASTA, myfile):
    for i in FASTA:
        print(i, file=myfile)
        print(FASTA[i], file=myfile)

path_input = os.path.join(os.getcwd(), "GENOMES")
path_output = os.path.join(os.getcwd(), 'FAA_changed_names')

if os.path.isdir(path_output):
    shutil.rmtree(path_output)
    os.mkdir(path_output)
else:
    os.mkdir(path_output)

print(path_output)

path_faa = os.path.join(path_input, "faa")
list_faa = os.listdir(path_faa)

path_gbff = os.path.join(path_input, "gbff")
list_gbff = os.listdir(path_gbff)

for j in list_gbff:
    path_j = os.path.join(path_gbff, j)
    assembly_name = re.compile('Assembly: .+')
    genome_name = re.compile('SOURCE      .+')
    with open(path_j, 'r') as file:
        gbff = file.read()
        G_name = re.search(genome_name, gbff).group().replace('SOURCE      ', '').replace('_', ' ')
        try:
            A_name = re.search(assembly_name, gbff).group().replace('Assembly: ', '')
        except:
            A_name = re.search(r'_[A-Z0-9]+v[0-9]+_', j).group().replace('_', '') # костыли на случай, если геном придётся скачивать из GENBANK, и в нём не окажется данных об аннотации

    for g in list_faa:
        if A_name in g:
            path_g = os.path.join(path_faa, g)
            with open(path_g, 'r') as file2:
                new_faa_list = fasta_for_analisys(file2.read(), G_name.replace(' ', '_'))
                new_file_name = G_name+'.faa'
                new_file_path = os.path.join(path_output, new_file_name)
                
            file3 = open(new_file_path, 'w')
            print_fasta(new_faa_list, file3)
            file3.close()
            print('the names of protein sequences have been rewritten for %s' % G_name)
