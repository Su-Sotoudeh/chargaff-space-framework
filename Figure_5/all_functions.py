import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

import gzip
from Bio import SeqIO
import pandas as pd

def opening_fasta_file(seq_path):
    fasta_path = seq_path
    if fasta_path.endswith('gz'):
        fasta_file = gzip.open(fasta_path,'rt')
    else:
        fasta_file = fasta_path
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(seq_record.seq)
        l  = len(sequence)/1000000
        #print('sequence id :',seq_record.id,' len chr :',l,'Mb')
        #print('opening_fasta_file finished')
        #records.append([str(seq_record.id),str(l),str(seq_record.seq)])

        #print(type(sequence))
    return sequence



def exon_and_intron_masking(sequence, gff_path,chr_id):

    #DNA  ATCG >> +1234 -5678
    #LTR  ATCG >> +LM@# -lm$%
    #SINE ATCG >> +OPJK -opjk
    #LINE ATCG >> +*YUI -+yui
    #Int  ATCG >> +QWER -qwer
    #Exe  ATCG >> +SDFH -sdfh
    #CDS  ATCG >> +ZXVB -zxvb
    #shre ATCG >> ATCG (shared between genes)
    #un   ATCG >> atcg



    chr_id = chr_id.split('.')[0]
    #if type(sequence) == list:
    #    #print(len(sequence))
    #    sequence = str(sequence[0])

    sequence = sequence.lower()  # ATCG all
    sequence_plus = sequence
    sequence_plus = list(sequence_plus) #ATCG all


#lableing genes***************************************************************************
    if gff_path.endswith('.gz'):
        gff = gzip.open(gff_path)
    else:
        gff = open(gff_path)
    df = pd.read_csv(gff,sep='\t',comment='#')
    df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    counter2 = 0
    ch_id_split = chr_id.split('.')[0]
    filtered_df = df[(df['seqid'].str.startswith(ch_id_split)) & (df['type'].str.startswith('gene'))]
    for annotation in filtered_df.to_numpy():

        counter2 += 1
        start = int(annotation[3])
        end = int(annotation[4])
        if end - start == 0:
            #print('len was 0')
            continue
        if annotation[6] == '+':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            gene_sequence = gene_sequence.replace('a','Q')
            gene_sequence = gene_sequence.replace('t','W')
            gene_sequence = gene_sequence.replace('c','E')
            gene_sequence = gene_sequence.replace('g','R')
            gene_sequence = gene_sequence.replace('q','A')
            gene_sequence = gene_sequence.replace('w','T')
            gene_sequence = gene_sequence.replace('e','C')
            gene_sequence = gene_sequence.replace('r','G')
            sequence_plus[(start - 1):end] = gene_sequence

        if annotation[6] == '-':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            gene_sequence = gene_sequence.replace('a','q')
            gene_sequence = gene_sequence.replace('t','w')
            gene_sequence = gene_sequence.replace('c','e')
            gene_sequence = gene_sequence.replace('g','r')
            gene_sequence = gene_sequence.replace('Q','A')
            gene_sequence = gene_sequence.replace('W','T')
            gene_sequence = gene_sequence.replace('E','C')
            gene_sequence = gene_sequence.replace('R','G')
            sequence_plus[(start - 1):end] = gene_sequence

    gff.close()
    gff.close()
    #print('gene pass',len(sequence_plus))

#lableing exons****************************************************************************

    if gff_path.endswith('.gz'):
        gff = gzip.open(gff_path)
    else:
        gff = open(gff_path)

    df = pd.read_csv(gff,sep='\t',comment='#')
    df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    counter2 = 0
    ch_id_split = chr_id.split('.')[0]
    filtered_df = df[(df['seqid'].str.startswith(ch_id_split)) & (df['type'].str.startswith('exon'))]
    for annotation in filtered_df.to_numpy():
            counter2 += 1
            start = int(annotation[3])
            end = int(annotation[4])
            if end - start == 0:
                #print('len was 0')
                continue




            if annotation[6] == '+':
                exon_sequence = sequence_plus[(start - 1):end]
                exon_sequence = ''.join(exon_sequence)
                exon_sequence = exon_sequence.replace('Q','S')
                exon_sequence = exon_sequence.replace('W','D')
                exon_sequence = exon_sequence.replace('E','F')
                exon_sequence = exon_sequence.replace('R','H')
                exon_sequence = exon_sequence.replace('q','A')
                exon_sequence = exon_sequence.replace('w','T')
                exon_sequence = exon_sequence.replace('e','C')
                exon_sequence = exon_sequence.replace('r','G')
                sequence_plus[(start - 1):end] = exon_sequence

            if annotation[6] == '-':
                exon_sequence = sequence_plus[(start - 1):end]
                exon_sequence = ''.join(exon_sequence)
                exon_sequence = exon_sequence.replace('q','s')
                exon_sequence = exon_sequence.replace('w','d')
                exon_sequence = exon_sequence.replace('e','f')
                exon_sequence = exon_sequence.replace('r','h')
                exon_sequence = exon_sequence.replace('Q','A')
                exon_sequence = exon_sequence.replace('W','T')
                exon_sequence = exon_sequence.replace('E','C')
                exon_sequence = exon_sequence.replace('R','G')
                sequence_plus[(start - 1):end] = exon_sequence
    gff.close()
    gff.close()
    #print('exon pass',len(sequence_plus))


#CDS ******************* CDS ********************* CDS
    if gff_path.endswith('.gz'):
        gff = gzip.open(gff_path)
    else:
        gff = open(gff_path)

    df = pd.read_csv(gff,sep='\t',comment='#')
    df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    counter2 = 0
    ch_id_split = chr_id.split('.')[0]
    filtered_df = df[(df['seqid'].str.startswith(ch_id_split)) & (df['type'].str.startswith('CDS'))]
    for annotation in filtered_df.to_numpy():
        counter2 += 1
        start = int(annotation[3])
        end = int(annotation[4])
        if end - start == 0:
            #print('len was 0')
            continue




        if annotation[6] == '+':
            CDS_sequence = sequence_plus[(start - 1):end]
            CDS_sequence = ''.join(CDS_sequence)
            CDS_sequence = CDS_sequence.replace('S','Z')
            CDS_sequence = CDS_sequence.replace('D','X')
            CDS_sequence = CDS_sequence.replace('F','V')
            CDS_sequence = CDS_sequence.replace('H','B')
            CDS_sequence = CDS_sequence.replace('z','A')
            CDS_sequence = CDS_sequence.replace('x','T')
            CDS_sequence = CDS_sequence.replace('v','C')
            CDS_sequence = CDS_sequence.replace('b','G')
            sequence_plus[(start - 1):end] = CDS_sequence

        if annotation[6] == '-':
            CDS_sequence = sequence_plus[(start - 1):end]
            CDS_sequence = ''.join(CDS_sequence)
            CDS_sequence = CDS_sequence.replace('s','z')
            CDS_sequence = CDS_sequence.replace('d','x')
            CDS_sequence = CDS_sequence.replace('f','v')
            CDS_sequence = CDS_sequence.replace('h','b')
            CDS_sequence = CDS_sequence.replace('Z','A')
            CDS_sequence = CDS_sequence.replace('X','T')
            CDS_sequence = CDS_sequence.replace('V','C')
            CDS_sequence = CDS_sequence.replace('B','G')
            sequence_plus[(start - 1):end] = CDS_sequence
    gff.close()
    gff.close()
    #print('CDS pass',len(sequence_plus))

    sequence_plus = ''.join(sequence_plus)
    #print('return joining',len(sequence_plus))
    return sequence_plus
def nucleotide_analysis_for_exon(sequence_plus):
    #un_sequenced_length = sequence_plus.count('N')  + sequence_plus.count('n')
    #sequenced_length = len(sequence_plus) - un_sequenced_length
    sequenced_length = len(sequence_plus)

    entron_P = round((sequence_plus.count('Q') + sequence_plus.count('W') + sequence_plus.count('E') + sequence_plus.count('R'))/sequenced_length*100,2)
    exon_p = round((sequence_plus.count('S') + sequence_plus.count('D') + sequence_plus.count('F') + sequence_plus.count('H'))/sequenced_length*100,2)
    CDS_P = round((sequence_plus.count('Z') + sequence_plus.count('X') + sequence_plus.count('V') + sequence_plus.count('B'))/sequenced_length*100,2)



    entron_C = round((sequence_plus.count('q') + sequence_plus.count('w') + sequence_plus.count('e') + sequence_plus.count('r'))/sequenced_length*100,2)
    exon_C = round((sequence_plus.count('s') + sequence_plus.count('d') + sequence_plus.count('f') + sequence_plus.count('h'))/sequenced_length*100,2)
    CDS_C = round((sequence_plus.count('z') + sequence_plus.count('x') + sequence_plus.count('v') + sequence_plus.count('b'))/sequenced_length*100,2)

    shered = round((sequence_plus.count('A') + sequence_plus.count('T') + sequence_plus.count('C') + sequence_plus.count('G'))/sequenced_length*100,2)


    #print('plus_information','entron',entron_P,'UTR',exon_p,'CDS',CDS_P)
    #print('comp_information','entron',entron_C,'UTR',exon_C,'CDS',CDS_C,)
    #print('masked_portion',(entron_P+exon_p+CDS_P+entron_C+exon_C+CDS_C),'+ Shered:',shered)


    A_entron_P = sequence_plus.count('Q')
    T_entron_P = sequence_plus.count('W')
    C_entron_P = sequence_plus.count('E')
    G_entron_P = sequence_plus.count('R')
    sum_entron_P = A_entron_P + T_entron_P + C_entron_P + G_entron_P
    Ap_entron_P = round(A_entron_P/sum_entron_P*100,3)
    Tp_entron_P = round(T_entron_P/sum_entron_P*100,3)
    Cp_entron_P = round(C_entron_P/sum_entron_P*100,3)
    Gp_entron_P = round(G_entron_P/sum_entron_P*100,3)

    A_entron_C = sequence_plus.count('q')
    T_entron_C = sequence_plus.count('w')
    C_entron_C = sequence_plus.count('e')
    G_entron_C = sequence_plus.count('r')
    sum_entron_C = A_entron_C + T_entron_C + C_entron_C + G_entron_C
    Ap_entron_C = round(A_entron_C/sum_entron_C*100,3)
    Tp_entron_C = round(T_entron_C/sum_entron_C*100,3)
    Cp_entron_C = round(C_entron_C/sum_entron_C*100,3)
    Gp_entron_C = round(G_entron_C/sum_entron_C*100,3)


#exon********************
    A_exon_P = sequence_plus.count('S')
    T_exon_P = sequence_plus.count('D')
    C_exon_P = sequence_plus.count('F')
    G_exon_P = sequence_plus.count('H')
    sum_exon_P = A_exon_P + T_exon_P + C_exon_P + G_exon_P
    Ap_exon_P = round(A_exon_P/sum_exon_P*100,3)
    Tp_exon_P = round(T_exon_P/sum_exon_P*100,3)
    Cp_exon_P = round(C_exon_P/sum_exon_P*100,3)
    Gp_exon_P = round(G_exon_P/sum_exon_P*100,3)

    A_exon_C = sequence_plus.count('s')
    T_exon_C = sequence_plus.count('d')
    C_exon_C = sequence_plus.count('f')
    G_exon_C = sequence_plus.count('h')
    sum_exon_C = A_exon_C + T_exon_C + C_exon_C + G_exon_C
    Ap_exon_C = round(A_exon_C/sum_exon_C*100,3)
    Tp_exon_C = round(T_exon_C/sum_exon_C*100,3)
    Cp_exon_C = round(C_exon_C/sum_exon_C*100,3)
    Gp_exon_C = round(G_exon_C/sum_exon_C*100,3)

    #CDS********************
    A_CDS_P = sequence_plus.count('Z')
    T_CDS_P = sequence_plus.count('X')
    C_CDS_P = sequence_plus.count('V')
    G_CDS_P = sequence_plus.count('B')
    sum_CDS_P = A_CDS_P + T_CDS_P + C_CDS_P + G_CDS_P
    Ap_CDS_P = round(A_CDS_P/sum_CDS_P*100,3)
    Tp_CDS_P = round(T_CDS_P/sum_CDS_P*100,3)
    Cp_CDS_P = round(C_CDS_P/sum_CDS_P*100,3)
    Gp_CDS_P = round(G_CDS_P/sum_CDS_P*100,3)

    A_CDS_C = sequence_plus.count('z')
    T_CDS_C = sequence_plus.count('x')
    C_CDS_C = sequence_plus.count('v')
    G_CDS_C = sequence_plus.count('b')
    sum_CDS_C = A_CDS_C + T_CDS_C + C_CDS_C + G_CDS_C
    Ap_CDS_C = round(A_CDS_C/sum_CDS_C*100,3)
    Tp_CDS_C = round(T_CDS_C/sum_CDS_C*100,3)
    Cp_CDS_C = round(C_CDS_C/sum_CDS_C*100,3)
    Gp_CDS_C = round(G_CDS_C/sum_CDS_C*100,3)

    import pandas as pd
    A = [Ap_entron_P,Ap_entron_C,Ap_exon_P,Ap_exon_C,Ap_CDS_P,Ap_CDS_C]
    T = [Tp_entron_P,Tp_entron_C,Tp_exon_P,Tp_exon_C,Tp_CDS_P,Tp_CDS_C]
    C = [Cp_entron_P,Cp_entron_C,Cp_exon_P,Cp_exon_C,Cp_CDS_P,Cp_CDS_C]
    G = [Gp_entron_P,Gp_entron_C,Gp_exon_P,Gp_exon_C,Gp_CDS_P,Gp_CDS_C]
    dic = {'A':A,'T':T,'C':C,'G':G}
    df = pd.DataFrame(dic)
    df.index = ['+entron','-entron','+exon','-exon','+CDS','-CDS']
    #print(df)

def TE_masking_after_exon(sequence_exon_masked,annotation_path,chr_id):

    #DNA  ATCG >> +1234 -5678
    #LTR  ATCG >> +LM@# -lm$%
    #SINE ATCG >> +OPJK -opjk
    #LINE ATCG >> +*YUI -+yui
    #Int  ATCG >> +QWER -qwer
    #Exe  ATCG >> +SDFH -sdfh
    #CDS  ATCG >> +ZXVB -zxvb
    #shre ATCG >> ATCG (shared between genes)
    #un   ATCG >> atcg


    sequence_plus = sequence_exon_masked
    sequence_plus = list(sequence_plus)


#lableing LINE***********************LINE***************************************LINE*************
    if annotation_path.endswith('.gz'):
        gff = gzip.open(annotation_path)
    else:
        gff = open(annotation_path)

    gff.readline()
    gff.readline()
    string_file = gff.read()
    string_file = string_file.replace(' *','*')
    from io import StringIO
    df = pd.read_csv(StringIO(string_file),sep='\\s+',comment='#')
    df.columns = ['SW score', '%div.', '%del.', '%ins.', 'seqid', 'begin', 'end', ' remain to end', 'orient','name','family','begin std','end','left','ID']
    df_filtered = df[(df['seqid'] ==chr_id)&(df['family'].str.startswith('LINE'))& (df['SW score']>=200)]





    for annotation in df_filtered.to_numpy():
        start = int(annotation[5])
        end = int(annotation[6])
        if end - start == 0:
            #print('len was 0')
            continue
        if annotation[8] == '+':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            replace_dict = str.maketrans("aQqSsZz", "*******")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("tWwDdXx", "YYYYYYY")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("cEeFfVv", "UUUUUUU")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("gRrHhBb", "IIIIIII")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            sequence_plus[(start - 1):end] = gene_sequence

        if annotation[8] == 'C':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            replace_dict = str.maketrans("aQqSsZz", "+++++++")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("tWwDdXx", "yyyyyyy")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("cEeFfVv", "uuuuuuu")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("gRrHhBb", "iiiiiii")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            sequence_plus[(start - 1):end] = gene_sequence
    gff.close()
    gff.close()

    #print('LINE pass',df_filtered.shape,'len chr',len(sequence_plus))
    #lableing SINE****************************SINE**************************************SINE*********
    if annotation_path.endswith('.gz'):
        gff = gzip.open(annotation_path)
    else:
        gff = open(annotation_path)

    gff.readline()
    gff.readline()
    string_file = gff.read()
    string_file = string_file.replace(' *','*')
    from io import StringIO
    df = pd.read_csv(StringIO(string_file),sep='\\s+',comment='#')
    df.columns = ['SW score', '%div.', '%del.', '%ins.', 'seqid', 'begin', 'end', ' remain to end', 'orient','name','family','begin std','end','left','ID']
    df_filtered = df[(df['seqid'] == chr_id)&(df['family'].str.startswith('SINE'))& (df['SW score']>=200)]
    for annotation in df_filtered.to_numpy():


        start = int(annotation[5])
        end = int(annotation[6])
        if end - start == 0:
            #print('len was 0')
            continue
        if annotation[8] == '+':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            replace_dict = str.maketrans("aQqSsZz", "OOOOOOO")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("tWwDdXx", "PPPPPPP")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("cEeFfVv", "JJJJJJJ")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("gRrHhBb", "KKKKKKK")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            sequence_plus[(start - 1):end] = gene_sequence

        if annotation[8] == 'C':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            replace_dict = str.maketrans("aQqSsZz", "ooooooo")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("tWwDdXx", "ppppppp")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("cEeFfVv", "jjjjjjj")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("gRrHhBb", "kkkkkkk")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            sequence_plus[(start - 1):end] = gene_sequence
    gff.close()
    gff.close()

    #print('SINE pass',df_filtered.shape,'len chr',len(sequence_plus))
    #lableing LTR************************LTR********************************************LTR*******
    if annotation_path.endswith('.gz'):
        gff = gzip.open(annotation_path)
    else:
        gff = open(annotation_path)

    gff.readline()
    gff.readline()
    string_file = gff.read()
    string_file = string_file.replace(' *','*')
    from io import StringIO
    df = pd.read_csv(StringIO(string_file),sep='\\s+',comment='#')
    df.columns = ['SW score', '%div.', '%del.', '%ins.', 'seqid', 'begin', 'end', ' remain to end', 'orient','name','family','begin std','end','left','ID']
    df_filtered = df[(df['seqid'] == chr_id)&(df['family'].str.startswith('LTR'))& (df['SW score']>=200)]
    for annotation in df_filtered.to_numpy():


        start = int(annotation[5])
        end = int(annotation[6])
        if end - start == 0:
            #print('len was 0')
            continue
        if annotation[8] == '+':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            replace_dict = str.maketrans("aQqSsZz", "LLLLLLL")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("tWwDdXx", "MMMMMMM")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("cEeFfVv", "@@@@@@@")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("gRrHhBb", "#######")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            sequence_plus[(start - 1):end] = gene_sequence

        if annotation[8] == 'C':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            replace_dict = str.maketrans("aQqSsZz", "lllllll")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("tWwDdXx", "mmmmmmm")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("cEeFfVv", "$$$$$$$")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("gRrHhBb", "%%%%%%%")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            sequence_plus[(start - 1):end] = gene_sequence
    gff.close()
    gff.close()
    #print('LTR pass',df_filtered.shape,'len chr',len(sequence_plus))
    #lableing DNA************************DNA********************************************DNA*******
    if annotation_path.endswith('.gz'):
        gff = gzip.open(annotation_path)
    else:
        gff = open(annotation_path)

    gff.readline()
    gff.readline()
    string_file = gff.read()
    string_file = string_file.replace(' *','*')
    from io import StringIO
    df = pd.read_csv(StringIO(string_file),sep='\\s+',comment='#')
    df.columns = ['SW score', '%div.', '%del.', '%ins.', 'seqid', 'begin', 'end', ' remain to end', 'orient','name','family','begin std','end','left','ID']
    df_filtered = df[(df['seqid'] == chr_id)&(df['family'].str.startswith('DNA'))& (df['SW score']>=200)]
    for annotation in df_filtered.to_numpy():
        start = int(annotation[5])
        end = int(annotation[6])
        if end - start == 0:
            #print('len was 0')
            continue
        if annotation[8] == '+':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            replace_dict = str.maketrans("aQqSsZz", "1111111")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("tWwDdXx", "2222222")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("cEeFfVv", "3333333")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("gRrHhBb", "4444444")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            sequence_plus[(start - 1):end] = gene_sequence

        if annotation[8] == 'C':
            gene_sequence = sequence_plus[(start - 1):end]
            gene_sequence = ''.join(gene_sequence)
            replace_dict = str.maketrans("aQqSsZz", "5555555")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("tWwDdXx", "6666666")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("cEeFfVv", "7777777")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            replace_dict = str.maketrans("gRrHhBb", "8888888")  # A, B, C → X
            gene_sequence = gene_sequence.translate(replace_dict)
            sequence_plus[(start - 1):end] = gene_sequence
    gff.close()
    gff.close()
    #print('DNA pass',df_filtered.shape,'len chr',len(sequence_plus))







    sequence_plus = ''.join(sequence_plus)
    #print(sequence_plus[:100])
    return sequence_plus
def nucleotide_analysis_for_exon_TE(sequence_plus):
    #un_sequenced_length = sequence_plus.count('N')  + sequence_plus.count('n')
    #sequenced_length = len(sequence_plus) - un_sequenced_length
    sequenced_length = len(sequence_plus)
    elements = []
    A_list = []
    T_list = []
    C_list = []
    G_list = []


    A_list_number = []
    T_list_number = []
    C_list_number = []
    G_list_number = []
    element_number_list = []
    # hi it is here *******************************************************************************
    #nN~~ have problem n=A N=T ~=C ~=G all nonsequence will be n as A
    list_letters = ['ZXVB','SDFH','QWER','*YUI','OPJK','LM@#','1234','zxvb','sdfh','qwer','+yui',
                    'opjk','lm$%','5678','ATCG','atcg','~~~~']#for ignoring ["nN~~"] converted  to ["~~~~"]






    for letters in list_letters:
        element = round((sequence_plus.count(letters[0]) + sequence_plus.count(letters[1]) + sequence_plus.count(letters[2]) + sequence_plus.count(letters[3]))/sequenced_length*100,3)
        element_in_number = sequence_plus.count(letters[0]) + sequence_plus.count(letters[1]) + sequence_plus.count(letters[2]) + sequence_plus.count(letters[3])

        A = sequence_plus.count(letters[0])
        T = sequence_plus.count(letters[1])
        C = sequence_plus.count(letters[2])
        G = sequence_plus.count(letters[3])
        sum_n = A + T + C + G
        try:
            Ap = round(A/sum_n*100,3)
            Tp = round(T/sum_n*100,3)
            Cp = round(C/sum_n*100,3)
            Gp = round(G/sum_n*100,3)
        except ZeroDivisionError:
            Ap = 0
            Tp = 0
            Cp = 0
            Gp = 0
        A_list.append(Ap)
        T_list.append(Tp)
        C_list.append(Cp)
        G_list.append(Gp)
        elements.append(element)
        
        A_list_number.append(A)
        T_list_number.append(T)
        C_list_number.append(C)
        G_list_number.append(G)
        element_number_list.append(element_in_number)





    A_list   =  A_list_number      +  [sum( A_list_number     )] + A_list
    T_list   =  T_list_number      +  [sum( T_list_number     )] + T_list
    C_list   =  C_list_number      +  [sum( C_list_number     )] + C_list
    G_list   =  G_list_number      +  [sum( G_list_number     )] + G_list
    elements = element_number_list +  [sum(element_number_list)] + elements





    import pandas as pd
    dic = {'A':A_list,'T':T_list,'C':C_list,'G':G_list,'whole':elements}
    df = pd.DataFrame(dic)
    df.index = ['CDS+','UTR+','Intron+','LINE+','SINE+','LTR+','DNA+','CDS-','UTR-','Intron-','LINE-','SINE-','LTR-','DNA-','Shared','Unannotated','unsequenced','sum','CDS+','UTR+','Intron+','LINE+','SINE+','LTR+','DNA+','CDS-','UTR-','Intron-','LINE-','SINE-','LTR-','DNA-','Shared','Unannotated','unsequenced']

         


    #print(df.to_string())
    #print('len chr',len(sequence_plus))
    '''number_list = []
    portion_list = []
    
    for i in "".join(list_letters):
        print(i,sequence_plus.count(i),sequence_plus.count(i)/len(sequence_plus*100))
        number_list.append(sequence_plus.count(i))
        portion_list.append(sequence_plus.count(i)/len(sequence_plus*100))
    print(sum(number_list),sum(portion_list))  '''
    '''plus = [0,2,4,6,8,10,12,14]
    minis = [1,3,5,7,9,11,13,14]
    P = 0
    M = 0
    for p,m in zip(plus,minis):
        P+=element_number_list[p]
        M+=element_number_list[m]

    print(P/M)'''

    '''entron_P = round((sequence_plus.count('Q') + sequence_plus.count('W') + sequence_plus.count('E') + sequence_plus.count('R'))/sequenced_length*100,2)
    exon_p = round((sequence_plus.count('S') + sequence_plus.count('D') + sequence_plus.count('F') + sequence_plus.count('H'))/sequenced_length*100,2)
    CDS_P = round((sequence_plus.count('Z') + sequence_plus.count('X') + sequence_plus.count('V') + sequence_plus.count('B'))/sequenced_length*100,2)



    entron_C = round((sequence_plus.count('q') + sequence_plus.count('w') + sequence_plus.count('e') + sequence_plus.count('r'))/sequenced_length*100,2)
    exon_C = round((sequence_plus.count('s') + sequence_plus.count('d') + sequence_plus.count('f') + sequence_plus.count('h'))/sequenced_length*100,2)
    CDS_C = round((sequence_plus.count('z') + sequence_plus.count('x') + sequence_plus.count('v') + sequence_plus.count('b'))/sequenced_length*100,2)

    shered = round((sequence_plus.count('A') + sequence_plus.count('T') + sequence_plus.count('C') + sequence_plus.count('G'))/sequenced_length*100,2)


    print('plus_information','entron',entron_P,'UTR',exon_p,'CDS',CDS_P)
    print('comp_information','entron',entron_C,'UTR',exon_C,'CDS',CDS_C,)
    print('masked_portion',(entron_P+exon_p+CDS_P+entron_C+exon_C+CDS_C),'+ Shered:',shered)


    A_entron_P = sequence_plus.count('Q')
    T_entron_P = sequence_plus.count('W')
    C_entron_P = sequence_plus.count('E')
    G_entron_P = sequence_plus.count('R')
    sum_entron_P = A_entron_P + T_entron_P + C_entron_P + G_entron_P
    Ap_entron_P = round(A_entron_P/sum_entron_P*100,3)
    Tp_entron_P = round(T_entron_P/sum_entron_P*100,3)
    Cp_entron_P = round(C_entron_P/sum_entron_P*100,3)
    Gp_entron_P = round(G_entron_P/sum_entron_P*100,3)

    A_entron_C = sequence_plus.count('q')
    T_entron_C = sequence_plus.count('w')
    C_entron_C = sequence_plus.count('e')
    G_entron_C = sequence_plus.count('r')
    sum_entron_C = A_entron_C + T_entron_C + C_entron_C + G_entron_C
    Ap_entron_C = round(A_entron_C/sum_entron_C*100,3)
    Tp_entron_C = round(T_entron_C/sum_entron_C*100,3)
    Cp_entron_C = round(C_entron_C/sum_entron_C*100,3)
    Gp_entron_C = round(G_entron_C/sum_entron_C*100,3)


    #exon********************
    A_exon_P = sequence_plus.count('S')
    T_exon_P = sequence_plus.count('D')
    C_exon_P = sequence_plus.count('F')
    G_exon_P = sequence_plus.count('H')
    sum_exon_P = A_exon_P + T_exon_P + C_exon_P + G_exon_P
    Ap_exon_P = round(A_exon_P/sum_exon_P*100,3)
    Tp_exon_P = round(T_exon_P/sum_exon_P*100,3)
    Cp_exon_P = round(C_exon_P/sum_exon_P*100,3)
    Gp_exon_P = round(G_exon_P/sum_exon_P*100,3)

    A_exon_C = sequence_plus.count('s')
    T_exon_C = sequence_plus.count('d')
    C_exon_C = sequence_plus.count('f')
    G_exon_C = sequence_plus.count('h')
    sum_exon_C = A_exon_C + T_exon_C + C_exon_C + G_exon_C
    Ap_exon_C = round(A_exon_C/sum_exon_C*100,3)
    Tp_exon_C = round(T_exon_C/sum_exon_C*100,3)
    Cp_exon_C = round(C_exon_C/sum_exon_C*100,3)
    Gp_exon_C = round(G_exon_C/sum_exon_C*100,3)

    #CDS********************
    A_CDS_P = sequence_plus.count('Z')
    T_CDS_P = sequence_plus.count('X')
    C_CDS_P = sequence_plus.count('V')
    G_CDS_P = sequence_plus.count('B')
    sum_CDS_P = A_CDS_P + T_CDS_P + C_CDS_P + G_CDS_P
    Ap_CDS_P = round(A_CDS_P/sum_CDS_P*100,3)
    Tp_CDS_P = round(T_CDS_P/sum_CDS_P*100,3)
    Cp_CDS_P = round(C_CDS_P/sum_CDS_P*100,3)
    Gp_CDS_P = round(G_CDS_P/sum_CDS_P*100,3)

    A_CDS_C = sequence_plus.count('z')
    T_CDS_C = sequence_plus.count('x')
    C_CDS_C = sequence_plus.count('v')
    G_CDS_C = sequence_plus.count('b')
    sum_CDS_C = A_CDS_C + T_CDS_C + C_CDS_C + G_CDS_C
    Ap_CDS_C = round(A_CDS_C/sum_CDS_C*100,3)
    Tp_CDS_C = round(T_CDS_C/sum_CDS_C*100,3)
    Cp_CDS_C = round(C_CDS_C/sum_CDS_C*100,3)
    Gp_CDS_C = round(G_CDS_C/sum_CDS_C*100,3)

    import pandas as pd
    A = [Ap_entron_P,Ap_entron_C,Ap_exon_P,Ap_exon_C,Ap_CDS_P,Ap_CDS_C]
    T = [Tp_entron_P,Tp_entron_C,Tp_exon_P,Tp_exon_C,Tp_CDS_P,Tp_CDS_C]
    C = [Cp_entron_P,Cp_entron_C,Cp_exon_P,Cp_exon_C,Cp_CDS_P,Cp_CDS_C]
    G = [Gp_entron_P,Gp_entron_C,Gp_exon_P,Gp_exon_C,Gp_CDS_P,Gp_CDS_C]
    dic = {'A':A,'T':T,'C':C,'G':G}
    df = pd.DataFrame(dic)
    df.index = ['+entron','-entron','+exon','-exon','+CDS','-CDS']
    print(df)'''

    return df


def analysis_datafram(df, path=None):
    print('hi')




    path = path or (config.MAMMALS_CHR_DIR / 'Human' / 'NC_000001.11.csv')
    import pandas as pd
    file = pd.read_csv(path)
    print(sum(file['whole']))
    print(file['whole'][15]/sum(file['whole'])*100)


















