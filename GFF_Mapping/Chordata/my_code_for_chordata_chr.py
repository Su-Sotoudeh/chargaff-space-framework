import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)

ON_SAVING_CSV = str(config.output_dir("GFF_Mapping", "Chordata") / 'on_saveing_time_resuts.csv')   # per-chromosome results (also used to resume)
import pandas as pd
import os
from multiprocessing import Pool
output_file = "output.csv"

for i in range(1, 6):  # example loop
    # Each iteration produces a new one-row DataFrame
    df = pd.DataFrame([[i * j for j in range(1, 11)]],
                      columns=[f"col{k}" for k in range(1, 11)])

    # Check if file exists and is non-empty
    write_header = not os.path.exists(output_file) or os.path.getsize(output_file) == 0

    # Append or write
    df.to_csv(output_file, mode='a', header=write_header, index=False)




# parallelized_version_of_your_script.py
import os
import gzip
import pandas as pd
import numpy as np
import time
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from Bio import SeqIO
import sys  # For potential memory estimation if needed




# ----------------- your original functions (kept intact, minimal fix) -----------------
def nucleotide_analysis_for_exon(sequence_plus):
    A_entron_P = sequence_plus.count('Q')
    T_entron_P = sequence_plus.count('W')
    C_entron_P = sequence_plus.count('E')
    G_entron_P = sequence_plus.count('R')
    sum_entron_P = A_entron_P + T_entron_P + C_entron_P + G_entron_P
    A_entron_C = sequence_plus.count('q')
    T_entron_C = sequence_plus.count('w')
    C_entron_C = sequence_plus.count('e')
    G_entron_C = sequence_plus.count('r')
    sum_entron_C = A_entron_C + T_entron_C + C_entron_C + G_entron_C
    A_exon_P = sequence_plus.count('S')
    T_exon_P = sequence_plus.count('D')
    C_exon_P = sequence_plus.count('F')
    G_exon_P = sequence_plus.count('H')
    sum_exon_P = A_exon_P + T_exon_P + C_exon_P + G_exon_P
    A_exon_C = sequence_plus.count('s')
    T_exon_C = sequence_plus.count('d')
    C_exon_C = sequence_plus.count('f')
    G_exon_C = sequence_plus.count('h')
    sum_exon_C = A_exon_C + T_exon_C + C_exon_C + G_exon_C

    A_CDS_P = sequence_plus.count('Z')
    T_CDS_P = sequence_plus.count('X')
    C_CDS_P = sequence_plus.count('V')
    G_CDS_P = sequence_plus.count('B')
    sum_CDS_P = A_CDS_P + T_CDS_P + C_CDS_P + G_CDS_P

    A_CDS_C = sequence_plus.count('z')
    T_CDS_C = sequence_plus.count('x')
    C_CDS_C = sequence_plus.count('v')
    G_CDS_C = sequence_plus.count('b')
    sum_CDS_C = A_CDS_C + T_CDS_C + C_CDS_C + G_CDS_C

    A_shared = sequence_plus.count('A')
    T_shared = sequence_plus.count('T')
    C_shared = sequence_plus.count('C')
    G_shared = sequence_plus.count('G')
    A_unannotated = sequence_plus.count('a')
    T_unannotated = sequence_plus.count('t')
    C_unannotated = sequence_plus.count('c')
    G_unannotated = sequence_plus.count('g')
    sum_unannotated = A_unannotated + T_unannotated + C_unannotated + G_unannotated

    # pandas import kept local to match your style (but imported at file top too)
    A = [A_entron_P, A_entron_C, A_exon_P, A_exon_C, A_CDS_P, A_CDS_C, A_shared, A_unannotated]
    T = [T_entron_P, T_entron_C, T_exon_P, T_exon_C, T_CDS_P, T_CDS_C, T_shared, T_unannotated]
    C = [C_entron_P, C_entron_C, C_exon_P, C_exon_C, C_CDS_P, C_CDS_C, C_shared, C_unannotated]
    G = [G_entron_P, G_entron_C, G_exon_P, G_exon_C, G_CDS_P, G_CDS_C, G_shared, G_unannotated]
    dic = {'A': A, 'T': T, 'C': C, 'G': G}
    df = pd.DataFrame(dic)
    df.index = ['+entron', '-entron', '+UTR', '-UTR', '+CDS', '-CDS', 'Shared', 'Unannotated']
    return df

def exon_and_intron_masking(sequence, df_all, chr_id):
    ch_id_split = chr_id.split('.')[0]
    sequence = sequence.lower()  # atcg all
    sequence_plus = bytearray(sequence.encode('ascii'))
    Inreon_P_counter = 0
    Inreon_C_counter = 0
    UTR_P_counter = 0
    UTR_C_counter = 0
    CDS_P_counter = 0
    CDS_C_counter = 0

    df_all.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    filtered_df_all = df_all[df_all['seqid'].str.startswith(ch_id_split)]

    # First pass: genes (introns)
    filtered_df = filtered_df_all[filtered_df_all['type'].str.startswith('gene')]
    for annotation in filtered_df.to_numpy():
        start = int(annotation[3])
        end = int(annotation[4])
        if end - start == 0:
            continue
        if annotation[6] == '+':
            Inreon_P_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('a'): sl[i] = ord('Q')
                elif ch == ord('t'): sl[i] = ord('W')
                elif ch == ord('c'): sl[i] = ord('E')
                elif ch == ord('g'): sl[i] = ord('R')
                elif ch == ord('q'): sl[i] = ord('A')
                elif ch == ord('w'): sl[i] = ord('T')
                elif ch == ord('e'): sl[i] = ord('C')
                elif ch == ord('r'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl
        if annotation[6] == '-':
            Inreon_C_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('a'): sl[i] = ord('q')
                elif ch == ord('t'): sl[i] = ord('w')
                elif ch == ord('c'): sl[i] = ord('e')
                elif ch == ord('g'): sl[i] = ord('r')
                elif ch == ord('Q'): sl[i] = ord('A')
                elif ch == ord('W'): sl[i] = ord('T')
                elif ch == ord('E'): sl[i] = ord('C')
                elif ch == ord('R'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl

    # Second pass: exons (UTRs)

    filtered_df = filtered_df_all[filtered_df_all['type'].str.startswith('exon')]
    for annotation in filtered_df.to_numpy():
        start = int(annotation[3])
        end = int(annotation[4])
        if end - start == 0:
            continue
        if annotation[6] == '+':
            UTR_P_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('Q'): sl[i] = ord('S')
                elif ch == ord('W'): sl[i] = ord('D')
                elif ch == ord('E'): sl[i] = ord('F')
                elif ch == ord('R'): sl[i] = ord('H')
                elif ch == ord('q'): sl[i] = ord('A')
                elif ch == ord('w'): sl[i] = ord('T')
                elif ch == ord('e'): sl[i] = ord('C')
                elif ch == ord('r'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl
        if annotation[6] == '-':
            UTR_C_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('q'): sl[i] = ord('s')
                elif ch == ord('w'): sl[i] = ord('d')
                elif ch == ord('e'): sl[i] = ord('f')
                elif ch == ord('r'): sl[i] = ord('h')
                elif ch == ord('Q'): sl[i] = ord('A')
                elif ch == ord('W'): sl[i] = ord('T')
                elif ch == ord('E'): sl[i] = ord('C')
                elif ch == ord('R'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl

    # Third pass: CDS
    filtered_df = filtered_df_all[filtered_df_all['type'].str.startswith('CDS')]
    for annotation in filtered_df.to_numpy():
        start = int(annotation[3])
        end = int(annotation[4])
        if end - start == 0:
            continue
        if annotation[6] == '+':
            CDS_P_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('S'): sl[i] = ord('Z')
                elif ch == ord('D'): sl[i] = ord('X')
                elif ch == ord('F'): sl[i] = ord('V')
                elif ch == ord('H'): sl[i] = ord('B')
                elif ch == ord('z'): sl[i] = ord('A')
                elif ch == ord('x'): sl[i] = ord('T')
                elif ch == ord('v'): sl[i] = ord('C')
                elif ch == ord('b'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl
        if annotation[6] == '-':
            CDS_C_counter += 1
            sl = sequence_plus[(start - 1):end]
            for i in range(len(sl)):
                ch = sl[i]
                if ch == ord('s'): sl[i] = ord('z')
                elif ch == ord('d'): sl[i] = ord('x')
                elif ch == ord('f'): sl[i] = ord('v')
                elif ch == ord('h'): sl[i] = ord('b')
                elif ch == ord('Z'): sl[i] = ord('A')
                elif ch == ord('X'): sl[i] = ord('T')
                elif ch == ord('V'): sl[i] = ord('C')
                elif ch == ord('B'): sl[i] = ord('G')
            sequence_plus[(start - 1):end] = sl

    sequence_plus = sequence_plus.decode('ascii')

    return sequence_plus, Inreon_P_counter, Inreon_C_counter, UTR_P_counter, UTR_C_counter, CDS_P_counter, CDS_C_counter

# ----------------- parallel worker: exact per-organism body kept -----------------
def process_one_chromosome(seq,organism,ch_id, gff_df):
    print('assined work')


    # Load GFF once per organism to save memory (avoid reloading multiple times per chromosome)



    my_table = None
    exon_mask_seq, Gene_P_counter, Gene_C_counter, Exon_P_counter, Exon_C_counter, CDS_P_counter, CDS_C_counter = exon_and_intron_masking(seq, gff_df, ch_id)
    if CDS_P_counter == 0 or CDS_C_counter == 0:
        print('CDS_P_counter == 0 or CDS_C_counter == 0')
        return None

    A_W = seq.count('A') + seq.count("a")
    T_W = seq.count('T') + seq.count("t")
    C_W = seq.count('C') + seq.count("c")
    G_W = seq.count("G") + seq.count("g")
    N_W = seq.count('N') + seq.count("n")
    Length = len(seq) - N_W
    organims_name = organism
    # protect division by zero (original assumes non-zero denominator). Keep behavior but safe:
    denom_W = (A_W + T_W + C_W)
    if denom_W == 0:
        GC_countent_list_W = 0
    else:
        GC_countent_list_W =(round((G_W + C_W) / denom_W * 100, 4))
    chr_id = ch_id
    len_seq = (round(Length / 1000000, 6))

    # avoid division by zero for T_W and C_W
    if T_W == 0:
        A2T = 0
    else:
        A2T = round((A_W) / (T_W), 4)
    if C_W == 0:
        G2C = 0
    else:
        G2C = (round((G_W) / (C_W), 4))

    L = len(exon_mask_seq)

    df = nucleotide_analysis_for_exon(exon_mask_seq)
    df.loc["Chr_Len"] = {"A": L, 'T': L, "C": L, "G": L, "whole": L}
    n_count = exon_mask_seq.count('n') + exon_mask_seq.count('N')
    df.loc["Nn_number"] = {"A": n_count, 'T': n_count, "C": n_count, "G": n_count, "whole": n_count}

    df.loc["organism"] = {"A": organism, 'T': organism, "C": organism, "G": organism, "whole": organism}
    df.loc["chr_id"] = {"A": ch_id, 'T': ch_id, "C": ch_id, "G": ch_id, "whole": ch_id}
    gene_elements = ['Int', 'UTR', 'CDS', 'Share', 'Unanot']
    counter = 0
    portions = []
    main_dict = {'Name': organims_name,
                 'Chr_ID': chr_id,
                 'Length(M)': len_seq,
                 'GC_W': GC_countent_list_W,
                 'A2T_W': A2T,
                 'G2C_W': G2C,}
    for e_number, element in enumerate(gene_elements):
        if element == "Share" or element == "Unanot":
            Ap, Tp, Cp, Gp = df.iloc[counter + 0]
            Ap, Tp, Cp, Gp = Ap / 2, Tp / 2, Cp / 2, Gp / 2
            Ac, Tc, Cc, Gc = df.iloc[counter + 0]
            Ac, Tc, Cc, Gc = Ac / 2, Tc / 2, Cc / 2, Gc / 2
            counter += 1
        else:
            Ap, Tp, Cp, Gp = df.iloc[counter + 0]
            Ac, Tc, Cc, Gc = df.iloc[counter + 1]
            counter += 2
        try:
            u = 1 / Ap + 1 / Tp + 1 / Cp + 1 / Gp
            u = 1 / Ac + 1 / Tc + 1 / Cc + 1 / Gc
            portions.append(sum([Ap, Tp, Cp, Gp, Ac, Tc, Cc, Gc]) / Length * 100)
            P_portion = round(Ap + Tp + Cp + Gp, 4)
            # protect division by zero inside these calculations; original code would raise;
            # we match behavior but return zeros when division by zero would occur
            if P_portion == 0:
                P_A_usage = 0
                P_T_usage = 0
                P_C_usage = 0
                P_G_usage =0
                P_Portion =0
                A2T_P = 0
                G2C_P = 0
                GC_countent_P = 0
            else:
                P_A_usage = (round(Ap / P_portion * 100, 4))
                P_T_usage = (round(Tp / P_portion * 100, 4))
                P_C_usage = (round(Cp / P_portion * 100, 4))
                P_G_usage = (round(Gp / P_portion * 100, 4))
                P_Portion = (round(P_portion / Length * 100, 4))
                # protect Tp==0
                if Tp == 0:
                    A2T_P= 0
                else:
                    A2T_P = (round(Ap / Tp, 4))
                if Cp == 0:
                    G2C_P = 0
                else:
                    G2C_P = (round(Gp / Cp, 4))

                GC_countent_P = (round((Cp + Gp) / P_portion * 100, 4))

            C_portion = round(Ac + Tc + Cc + Gc, 4)
            if C_portion == 0:
                C_A_usage = 0
                C_T_usage = 0
                C_C_usage = 0
                C_G_usage = 0
                C_portion = 0
                A2T_C = 0
                G2C_C = 0
                GC_countent_C = 0
            else:
                C_A_usage=(round(Ac / C_portion * 100, 4))
                C_T_usage=(round(Tc / C_portion * 100, 4))
                C_C_usage=(round(Cc / C_portion * 100, 4))
                C_G_usage=(round(Gc / C_portion * 100, 4))
                C_portion=(round(C_portion / Length * 100, 4))
                if Tc == 0:
                    A2T_C=0
                else:
                    A2T_C=(round(Ac / Tc, 4))
                if Cc == 0:
                    G2C_C=0
                else:
                    G2C_C = (round(Gc / Cc, 4))
                C_element_number = (CDS_C_counter)
                GC_countent_C = (round((Cc + Gc) / C_portion * 100, 4))

            A2T_anot = (round((Ap + Ac) / (Tp + Tc) if (Tp + Tc) != 0 else 0, 4))
            G2C_anot = (round((Gp + Gc) / (Cp + Cc) if (Cp + Cc) != 0 else 0, 4))

            P2C = (round(P_portion / C_portion, 4) if C_portion != 0 else 0)

            PandC = [P_portion, C_portion]
            PandC.sort()
            if PandC[1] == 0:
                smaller_portion_2_Bigger = 0
            else:
                smaller_portion_2_Bigger = (round(PandC[0] / PandC[1], 4))
            dif2sum = (round(abs(PandC[0] - PandC[1]) / ((PandC[0] + PandC[1]) * 2) * 100, 4) if (PandC[0] + PandC[1]) != 0 else 0)

        except ZeroDivisionError:
            portions.append(0)
            P_A_usage = 0
            P_T_usage = 0
            P_C_usage = 0
            P_G_usage = 0
            P_Portion = 0
            A2T_P  = 0
            G2C_P  = 0

            GC_countent_P = 0

            P_portion = 0
            C_A_usage = 0
            C_T_usage = 0
            C_C_usage = 0
            C_G_usage = 0
            C_portion = 0
            A2T_C = 0
            G2C_C = 0
            GC_countent_C = 0
            A2T_anot = 0
            G2C_anot = 0
            P2C = 0
            PandC = [P_portion, C_portion]
            PandC.sort()
            smaller_portion_2_Bigger = 0
            dif2sum = 0

        vaiavle_dic = {
            f'{element}_GC_P': GC_countent_P,
            f'{element}_GC_C': GC_countent_C,
            f'{element}_P%': P_Portion,
            f'{element}_C%': C_portion,
            f'{element}_P2C': P2C,
            f'{element}_S2B': smaller_portion_2_Bigger,
            f'{element}_dif2sum': dif2sum,
            f'{element}_A2T_ano': A2T_anot,
            f'{element}_G2C_ano': G2C_anot,
            f'{element}_P_A_ano': P_A_usage,
            f'{element}_P_T_ano': P_T_usage,
            f'{element}_P_C_ano': P_C_usage,
            f'{element}_P_G_ano': P_G_usage,
            f'{element}_P_A2T':   A2T_P,
            f'{element}_P_G2C': G2C_P,
            f'{element}_C_A_ano': C_A_usage,
            f'{element}_C_T_ano': C_T_usage,
            f'{element}_C_C_ano': C_C_usage,
            f'{element}_C_G_ano': C_G_usage,
            f'{element}_A2T_C': A2T_C,
            f'{element}_G2C_C': G2C_C
        }

        main_dict = main_dict | vaiavle_dic


    All_annotation = (round(sum(portions[:4]), 4))
    Int_annotaion = (round(portions[0], 4))
    UTR_annotation = (round(portions[1], 4))
    CDS_annotation = (round(portions[2], 4))
    Shared_annotation = (round(portions[3], 4))
    Un_annotation = (round(portions[4], 4))
    Sum_E = (round(sum(portions), 4))
    try:
        my_table = pd.DataFrame([main_dict])

    except ValueError:
        for i in main_dict:
            print('hi', i, main_dict[i])
    my_table.insert(4, 'Wh_ano', All_annotation)
    my_table.insert(5, 'Int_ano', Int_annotaion)
    my_table.insert(6, 'UTR_ano', UTR_annotation)
    my_table.insert(7, 'CDS_ano', CDS_annotation)
    my_table.insert(8, 'share_ano', Shared_annotation)
    my_table.insert(9, 'Un_ano', Un_annotation)
    my_table.insert(10, 'sum', Sum_E)
    my_table['Gene_P_count'] = Gene_P_counter
    my_table['Gene_C_count'] = Gene_C_counter
    my_table['Exon_P_count'] = Exon_P_counter
    my_table['Exon_C_count'] = Exon_C_counter
    my_table['CDS_P_count'] = CDS_P_counter
    my_table['CDS_C_count'] = CDS_C_counter

    column_width = 15
    header = " ".join(col.ljust(column_width) for col in my_table.columns)
    last_row_values = " ".join(str(val).ljust(column_width) for val in my_table.iloc[-1].values)
    '''if analysed_chr_counter == 1:
        print(header)
        print("-" * len(header)) '''
    print(last_row_values)
    #here it should be saved----------------------------------------------------------------------
    #write_header = not os.path.exists('on_saveing_time_resuts_new.csv') or os.path.getsize('on_saveing_time_resuts_new.csv') == 0
    #my_table.to_csv(os.path.join('on_saveing_time_resuts_new.csv'), mode='a', header=write_header, index=False)
    return  my_table





# ----------------- parallel controller (keeps your outer behavior but parallel) -----------------
def prepar_chromosomes(all_genomes_folder_path,organ_group):
    print('inalysis_started')
    start_time = time.time()
    current_time = datetime.now()
    start_run = current_time.strftime("%Y-%m-%d_%H-%M-%S")
    print(start_run)
    organism_list = os.listdir(all_genomes_folder_path)
    organism_number = len(organism_list)
    counter = 0
    import pandas as pd
    # resume: reuse results of an earlier run if they exist, otherwise start from scratch
    saved_file = pd.read_csv(ON_SAVING_CSV) if os.path.exists(ON_SAVING_CSV) else pd.DataFrame(columns=['Name', 'Chr_ID'])
    the_last_organ = saved_file.tail(8)['Name'].tolist()
    print('the_last_organ',the_last_organ)
    saved_chr_ids = saved_file['Chr_ID'].tolist()
    saved_names = saved_file['Name'].tolist()
    for i in the_last_organ:
        saved_names.remove(i)
    #print(len(organism_list),organism_list)

    for organism in organism_list:
        counter += 1
        print(organism_number,counter,organism)
        if counter < 625 :
            continue

        if organism in saved_names  :
            print(organism,'before analysed')
            continue
        if organism.endswith(".json") or organism.endswith(".jsonl"):
            continue
        files = os.listdir(all_genomes_folder_path + '/' + organism)
        fna = [i for i in files if i.endswith('.fna.gz') and not i.startswith('cds')]
        gff = [k for k in files if k.endswith('.gff.gz')]
        if len(gff) == 0 or len(fna) == 0:
            print('len(gff) == 0 or len(fna) == 0')
            continue
        fna.sort()
        gff.sort()
        try:
            fasta_file = gzip.open(all_genomes_folder_path + '/' + organism + '/' + fna[0], 'rt')
        except (IOError, EOFError, IndexError):
            continue
        # Load GFF once per organism to save memory (avoid reloading multiple times per chromosome)
        gff_path_main = all_genomes_folder_path + '/' + organism + '/' + gff[0]
        try:
            if gff_path_main.endswith('.gz'):
                with gzip.open(gff_path_main, 'rt') as gff_file:
                    df_all = pd.read_csv(gff_file, sep='\t', comment='#')
            else:
                with open(gff_path_main, 'r') as gff_file:
                    df_all = pd.read_csv(gff_file, sep='\t', comment='#')
        except (IOError, EOFError, IndexError):
            continue

        #here is a test------------------------------------------------------------------------------------------------------
        chromosomes_pool = []
        chr_id_pool = []
        organ_name_pool = []
        gff_path_pool = []
        try:
            for seq_record in SeqIO.parse(fasta_file, "fasta"):


                ch_id = seq_record.id
                if ch_id in saved_chr_ids:
                    print(ch_id," this chromosome already analysed")
                    continue
                if ch_id in saved_chr_ids:
                    print('chr_exit_and analysed')
                    continue
                if 'mitochondrion' in seq_record.description:
                    continue
                seq = str(seq_record.seq)
                l = len(seq) / 1000000
                if l < 0.05 or l >200:
                    continue
                if "chromosome" not in seq_record.description:
                    continue
                if ch_id.startswith('NW'):
                    continue
                chromosomes_pool.append(seq)
                chr_id_pool.append(ch_id)
                organ_name_pool.append(organism)
                gff_path_pool.append(df_all)

                if len(chromosomes_pool) == 8:
                    tasks = [(seq, organ,id ,gff_path) for seq,organ,id,gff_path in zip(chromosomes_pool,organ_name_pool,chr_id_pool,gff_path_pool)]
                    with Pool(processes=8) as pool:
                        results = pool.starmap(process_one_chromosome, tasks)
                        df_final = pd.concat(results, ignore_index=True)
                        write_header = not os.path.exists(ON_SAVING_CSV) or os.path.getsize(ON_SAVING_CSV) == 0
                        df_final.to_csv(ON_SAVING_CSV, mode='a', header=write_header, index=False)
                    chromosomes_pool  = []
                    chr_id_pool       = []
                    organ_name_pool   = []
                    gff_path_pool     = []
        except EOFError:
            print('why')
            continue
    end_time = time.time()
    print('procces time', round((end_time - start_time) / 60, 2), 'minit')
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    print('start running', start_run)
    print('end running  ', timestamp)
    print('procces finished')

genomes_path = str(config.require(config.ncbi_data_dir("Chordata"), "Unzip the NCBI genome package of this group into Fasta_GFF_inputs/<Group>/ (see Fasta_GFF_inputs/README.md)."))   # Fasta_GFF_inputs/Chordata/ncbi_dataset/data

prepar_chromosomes(genomes_path,'Chordata')





