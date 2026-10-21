import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)
import os
import pandas as pd

one_Trabite = str(config.require(config.ncbi_data_dir("Fungi"), "Unzip the NCBI genome package of this group into Fasta_GFF_inputs/<Group>/ (see Fasta_GFF_inputs/README.md).")) + "/"   # Fasta_GFF_inputs/Fungi/ncbi_dataset/data/
organism_names = os.listdir(one_Trabite)

refrence_list = pd.read_csv(config.require(config.FUNGI_REFERENCE_LIST, "Copy fungi_ref_annotated_2025_10_14.tsv into Dataframes/."), sep='\t')
refrence_list = list(refrence_list['Assembly Accession'])


def gene_element_annotation_in_plus_and_minis(all_genomes_folder_path,refrence_list):
    print('inalysis_started')
    # listing file
    import os
    import gzip
    import pandas as pd
    import numpy as np
    import time
    start_time = time.time()
    from datetime import datetime
    current_time = datetime.now()
    # Format it as a string for naming files (e.g., "2024-12-15_14-30-45")
    start_run = current_time.strftime("%Y-%m-%d_%H-%M-%S")
    print(start_run)
    #Intron
    A2T_list =                      []
    G2C_list =                      []
    P2C_number_list =               [[],[],[],[],[]]
    GC_countent_list_W =            []
    GC_countent_list_P =            [[],[],[],[],[]]
    GC_countent_list_C =            [[],[],[],[],[]]
    len_seq_list =                  []
    plus_and_minis_strand_portion = [[],[],[],[],[]]
    organims_name =                 []
    A2T_anot_list =                 [[],[],[],[],[]]
    G2C_anot_list =                 [[],[],[],[],[]]
    P_element_number_list =         [[],[],[],[],[]]
    C_element_number_list =         [[],[],[],[],[]]

    P_A_usage_list =                [[],[],[],[],[]]
    P_T_usage_list =                [[],[],[],[],[]]
    P_C_usage_list =                [[],[],[],[],[]]
    P_G_usage_list =                [[],[],[],[],[]]

    C_A_usage_list =                [[],[],[],[],[]]
    C_T_usage_list =                [[],[],[],[],[]]
    C_C_usage_list =                [[],[],[],[],[]]
    C_G_usage_list =                [[],[],[],[],[]]
    smaller_portion_2_Bigger =      [[],[],[],[],[]]
    dif2sum =                       [[],[],[],[],[]]
    P_Portion_list =                [[],[],[],[],[]]
    C_portion_list =                [[],[],[],[],[]]
    P2C_list =                      [[],[],[],[],[]]
    P2C_list_W =                    [[],[],[],[],[]]

    A2T_P_list =                    [[],[],[],[],[]]
    G2C_P_list =                    [[],[],[],[],[]]
    A2T_C_list =                    [[],[],[],[],[]]
    G2C_C_list =                    [[],[],[],[],[]]

    chr_id_list =                   []
    All_annotation =                []
    Int_annotaion =                 []
    UTR_annotation =                []
    CDS_annotation =                []
    Shared_annotation =             []
    Un_annotation =                 []
    Sum_E =                         []

    Gene_P_counter_list  =               []
    Gene_C_counter_list  =               []
    Exon_P_counter_list  =               []
    Exon_C_counter_list  =               []
    CDS_P_counter_list  =               []
    CDS_C_counter_list  =               []


    analysed_chr_counter = 0
    lablename = 'refrence'
    organ_group = "fungis"


    #print(organisms)
    organism_list = os.listdir(all_genomes_folder_path)
    #for organism in organisms:
    organism_counter = 0
    #print(organism_list)


    for organism in organism_list: # 1 time use not in loop
        #if organism_counter  == 5:
        #    break
        if organism not in refrence_list:
            continue

        if organism.endswith(".json") or organism.endswith(".jsonl"):
           continue
        organism_counter+=1
        #print(organism,' analysis')
        files = os.listdir(all_genomes_folder_path+'/'+organism)
        #print(files)

        fna = [i for i in files if i.endswith('.fna.gz') and not i.startswith('cds')]
        #out = [j for j in files if j.endswith('.out')]
        gff = [k for k in files if k.endswith('.gff.gz')]
        csv = [p for p in files if p.endswith('.csv')]
        fna.sort()
        #out.sort()
        gff.sort()
        from all_functions import  exon_and_intron_masking
        from all_functions import opening_fasta_file
        from all_functions import nucleotide_analysis_for_exon
        from Bio import SeqIO


        try:
            fasta_file = gzip.open(all_genomes_folder_path+'/'+organism+'/'+fna[0],'rt')
        except IOError or EOFError:
            continue
        try:
            analysed_chr_counter = 0
            for seq_record in SeqIO.parse(fasta_file, "fasta"):
                if 'mitochondrion' in seq_record.description:
                    continue
                    #print(seq_record.description)
                analysed_chr_counter += 1
                seq = str(seq_record.seq)
                l  = len(seq)/1000000
                if l < 0.05:
                    continue
                if "chromosome" not in seq_record.description:
                    continue
                ch_id = seq_record.id
                chr_number = (seq_record.description.split())[0]
                #print('sequence id :',ch_id,' len chr :',l,'Mb')
                #print('opening_fasta_file finished')





                #print('A',A_W,'T',T_W,'C',C_W,'G',G_W,'N',N_W,'Lenght',Length)
                exon_mask_seq,Gene_P_counter,Gene_C_counter,Exon_P_counter,Exon_C_counter ,CDS_P_counter,CDS_C_counter  = exon_and_intron_masking(seq,all_genomes_folder_path+'/'+organism+'/'+gff[0],ch_id)
                if CDS_P_counter == 0 or CDS_C_counter == 0:
                    continue


                Gene_P_counter_list.append(Gene_P_counter)
                Gene_C_counter_list.append(Gene_C_counter)
                Exon_P_counter_list.append(Exon_P_counter)
                Exon_C_counter_list.append(Exon_C_counter)
                CDS_P_counter_list.append(CDS_P_counter)
                CDS_C_counter_list.append(CDS_C_counter)

                A_W =  seq.count('A')+seq.count("a")
                T_W =  seq.count('T')+seq.count("t")
                C_W =  seq.count('C')+seq.count("c")
                G_W =  seq.count("G")+seq.count("g")
                N_W =  seq.count('N')+seq.count("n")
                Length = len(seq) - N_W
                organims_name.append(organism)
                GC_countent_list_W.append(round((G_W+C_W)/(A_W+T_W+C_W)*100,4))
                chr_id_list.append(ch_id)
                len_seq_list.append(round(Length/1000000,6))
                A2T_list.append(round((A_W)/(T_W),4))
                G2C_list.append(round((G_W)/(C_W),4))

                #print('exon_and_intron_masking finished')

                L = len(exon_mask_seq)

                df = nucleotide_analysis_for_exon(exon_mask_seq)
                #N = 1-(sum(df.iloc[7]))/L*100
                #print(N)
                df.loc["Chr_Len"] = {"A":L,'T':L,"C":L,"G":L,"whole":L}
                n_count = exon_mask_seq.count('n')+exon_mask_seq.count('N')
                df.loc["Nn_number"] = {"A":n_count,'T':n_count,"C":n_count,"G":n_count,"whole":n_count}

                df.loc["organism"] = {"A":organism,'T':organism,"C":organism,"G":organism,"whole":organism}
                df.loc["chr_id"] = {"A":ch_id,'T':ch_id,"C":ch_id,"G":ch_id,"whole":ch_id}
                #print(df)


                gene_elements = ['Int','UTR','CDS','Share','Unanot']
                counter = 0
                portions = []
                main_dict={'Name':organims_name,
                         'Chr_ID':chr_id_list,
                         'Length(M)':len_seq_list,
                         'GC_W':GC_countent_list_W,
                         'A2T_W':A2T_list,
                         'G2C_W':G2C_list,}
                for e_number,element in enumerate(gene_elements):
                    #print(element)
                    if element == "Share" or element == "Unanot":

                        Ap,Tp,Cp,Gp = df.iloc[counter+0]
                        Ap,Tp,Cp,Gp = Ap/2,Tp/2,Cp/2,Gp/2
                        Ac,Tc,Cc,Gc = df.iloc[counter+0]
                        Ac,Tc,Cc,Gc = Ac/2,Tc/2,Cc/2,Gc/2
                        counter += 1
                        #print(element,'+',Ap,Tp,Cp,Gp)
                        #print(element,'-',Ac,Tc,Cc,Gc)
                    else:

                        Ap,Tp,Cp,Gp = df.iloc[counter+0]
                        Ac,Tc,Cc,Gc = df.iloc[counter+1]
                        #print(element,'+',Ap,Tp,Cp,Gp)
                        #print(element,'-',Ac,Tc,Cc,Gc)

                        counter += 2

                    try :
                        u = 1/Ap + 1/Tp + 1/Cp + 1/Gp
                        u = 1/Ac + 1/Tc + 1/Cc + 1/Gc
                        portions.append(sum([Ap,Tp,Cp,Gp,Ac,Tc,Cc,Gc])/Length*100)

                        #ppppppppppppppppppppppppppppppppp
                        P_portion = round(Ap+Tp+Cp+Gp,4)
                        P_A_usage_list[e_number].append(round(Ap/P_portion*100,4))
                        P_T_usage_list[e_number].append(round(Tp/P_portion*100,4))
                        P_C_usage_list[e_number].append(round(Cp/P_portion*100,4))
                        P_G_usage_list[e_number].append(round(Gp/P_portion*100,4))
                        P_Portion_list[e_number].append(round(P_portion/Length*100,4))
                        A2T_P_list[e_number].append(round(Ap/Tp,4))
                        G2C_P_list[e_number].append(round(Gp/Cp,4))
                        P_element_number_list[e_number].append(CDS_P_counter)
                        GC_countent_list_P[e_number].append(round((Cp+Gp)/P_portion*100,4))

                        #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                        C_portion = round(Ac+Tc+Cc+Gc,4)
                        C_A_usage_list[e_number].append(round(Ac/C_portion*100,4))
                        C_T_usage_list[e_number].append(round(Tc/C_portion*100,4))
                        C_C_usage_list[e_number].append(round(Cc/C_portion*100,4))
                        C_G_usage_list[e_number].append(round(Gc/C_portion*100,4))
                        C_portion_list[e_number].append(round(P_portion/Length*100,4))
                        A2T_C_list[e_number].append(round(Ac/Tc,4))
                        G2C_C_list[e_number].append(round(Gc/Cc,4))
                        C_element_number_list[e_number].append(CDS_C_counter)
                        GC_countent_list_C[e_number].append(round((Cc+Gc)/C_portion*100,4))

                        #annotat
                        A2T_anot_list[e_number].append(round((Ap+Ac)/(Tp+Tc),4))
                        G2C_anot_list[e_number].append(round((Gp+Gc)/(Cp+Cc),4))


                        P2C_list[e_number].append(round(P_portion/C_portion,4))
                        P2C_number_list[e_number].append(round(CDS_P_counter/CDS_C_counter,4))
                        PandC = [P_portion,C_portion]
                        PandC.sort()
                        smaller_portion_2_Bigger[e_number].append(round(PandC[0]/PandC[1],4))
                        dif2sum[e_number].append(round(abs(PandC[0]-PandC[1])/((PandC[0]+PandC[1])*2)*100,4))
                    except ZeroDivisionError:

                        #print(element,'*****',Ap,Tp,Cp,Gp)

                        portions.append(0)
                        #ppppppppppppppppppppppppppppppppp
                        P_portion = round(Ap+Tp+Cp+Gp,4)
                        P_A_usage_list[e_number].append(0)
                        P_T_usage_list[e_number].append(0)
                        P_C_usage_list[e_number].append(0)
                        P_G_usage_list[e_number].append(0)
                        P_Portion_list[e_number].append(0)
                        A2T_P_list[e_number].append(0)
                        G2C_P_list[e_number].append(0)
                        P_element_number_list[e_number].append(0)
                        GC_countent_list_P[e_number].append(0)

                        #CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                        C_portion = round(Ac+Tc+Cc+Gc,4)
                        C_A_usage_list[e_number].append(0)
                        C_T_usage_list[e_number].append(0)
                        C_C_usage_list[e_number].append(0)
                        C_G_usage_list[e_number].append(0)
                        C_portion_list[e_number].append(0)
                        A2T_C_list[e_number].append(0)
                        G2C_C_list[e_number].append(0)
                        C_element_number_list[e_number].append(0)
                        GC_countent_list_C[e_number].append(0)

                        #annotat
                        A2T_anot_list[e_number].append(0)
                        G2C_anot_list[e_number].append(0)

                        P2C_list[e_number].append(0)
                        P2C_number_list[e_number].append(0)
                        PandC = [P_portion,C_portion]
                        PandC.sort()
                        smaller_portion_2_Bigger[e_number].append(0)
                        dif2sum[e_number].append(0)



                    vaiavle_dic={
                                 f'{element}_GC_P':GC_countent_list_P[e_number],
                                 f'{element}_GC_C':GC_countent_list_C[e_number],
                                 f'{element}_P%':P_Portion_list[e_number],
                                 f'{element}_C%':C_portion_list[e_number],
                                 f'{element}_P2C':P2C_list[e_number],
                                 f'{element}_S2B':smaller_portion_2_Bigger[e_number],
                                 f'{element}_dif2sum':dif2sum[e_number],
                                 f'{element}_A2T_ano':A2T_anot_list[e_number],
                                 f'{element}_G2C_ano':G2C_anot_list[e_number],
                                 f'{element}_P_E_num':P_element_number_list[e_number],
                                 f'{element}_C_E_num':C_element_number_list[e_number],
                                 f'{element}_P_A_ano':P_A_usage_list[e_number],
                                 f'{element}_P_T_ano':P_T_usage_list[e_number],
                                 f'{element}_P_C_ano':P_C_usage_list[e_number],
                                 f'{element}_P_G_ano':P_G_usage_list[e_number],
                                 f'{element}_P_A2T'  :A2T_P_list[e_number],
                                 f'{element}_P_G2C'  :G2C_P_list[e_number],
                                 f'{element}_C_A_ano':C_A_usage_list[e_number],
                                 f'{element}_C_T_ano':C_T_usage_list[e_number],
                                 f'{element}_C_C_ano':C_C_usage_list[e_number],
                                 f'{element}_C_G_ano':C_G_usage_list[e_number],
                                 f'{element}_A2T_C':A2T_C_list[e_number],
                                 f'{element}_G2C_C':G2C_C_list[e_number]}

                    main_dict = main_dict | vaiavle_dic
                All_annotation.append(round(sum(portions[:4]),4))
                Int_annotaion.append(round(portions[0],4))
                UTR_annotation.append(round(portions[1],4))
                CDS_annotation.append(round(portions[2],4))
                Shared_annotation.append(round(portions[3],4))
                Un_annotation.append(round(portions[4],4))
                Sum_E.append(round(sum(portions),4))
                #print(sum(portions[:4]))
                #main_dict.inser[f'All_anno'] = All_annotation
                try:
                    my_table = pd.DataFrame(main_dict)
                except ValueError:
                    for i in main_dict:
                        print('hi',i,len(main_dict[i]))
                my_table.insert(4,'Wh_ano',All_annotation)
                my_table.insert(5,'Int_ano',Int_annotaion)
                my_table.insert(6,'UTR_ano',UTR_annotation)
                my_table.insert(7,'CDS_ano',CDS_annotation)
                my_table.insert(8,'share_ano',Shared_annotation)
                my_table.insert(9,'Un_ano',Un_annotation)
                my_table.insert(10,'sum',Sum_E)
                my_table['Gene_P_count'] = Gene_P_counter_list
                my_table['Gene_C_count'] = Gene_C_counter_list
                my_table['Exon_P_count'] = Exon_P_counter_list
                my_table['Exon_C_count'] = Exon_C_counter_list
                my_table['CDS_P_count']  = CDS_P_counter_list
                my_table['CDS_C_count']  = CDS_C_counter_list



                # Print the last row with aligned columns
                column_width = 15# Adjust width as needed
                header = " ".join(col.ljust(column_width) for col in my_table.columns)
                last_row_values = " ".join(str(val).ljust(column_width) for val in my_table.iloc[-1].values)

                # Print the header once
                #if len(my_table) == 1:  # If it's the first row being added
                if analysed_chr_counter == 1:
                    print(header)
                    print("-" * len(header))  # Add a separator line for clarity

                # Print the last row in aligned format

                print(organism_counter,analysed_chr_counter,last_row_values)
                #print(sum(portions[:4]),portions)

        except EOFError:
            continue
    import time
    end_time = time.time()

    print('procces time',round((end_time-start_time)/60,2),'minit')

    current_time = datetime.now()

    # Format it as a string for naming files (e.g., "2024-12-15_14-30-45")
    timestamp = current_time.strftime("%Y-%m-%d_%H-%M-%S")

    results_folder = os.path.join(str(config.output_dir("GFF_Mapping", organ_group)), f"{organ_group}{lablename}_Int_UTR_CDS_{timestamp}")
    os.makedirs(results_folder, exist_ok=True)  # Creates the folder if it doesn't exist

    # Define the output file path
    output_filename = os.path.join(results_folder, f"results_{timestamp}.xlsx")

    # Save the DataFrame to the timestamped Excel file
    my_table.to_excel(output_filename, index=False)

    print(f"Results saved to: {output_filename}")
    print('start running',start_run)
    print('end running  ',timestamp)
    print('procces finished')



gene_element_annotation_in_plus_and_minis(one_Trabite,refrence_list)













