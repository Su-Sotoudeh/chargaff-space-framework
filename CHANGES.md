# What changed (paths and folder names only; no analysis code was changed)

config.py
- Masked_Fastas -> **Masked_Fasta**; Dataframes/Mammals -> **Dataframes/Mammals_CSV**
- new: MAMMALS_WORK_DIR (= Fasta_GFF_inputs/Mammals), REPEATMASKER_BIN (env CSPR_REPEATMASKER), REPEATMASKER_THREADS (env CSPR_RM_THREADS), figure_input(n, name)
- removed the unused VIRUS_CSV, VIRUS_REF_CSV, PLASMID_CSV

all_functions.py moved from Figure_5/ to the repository root (all scripts find it through config.py).
Figure_3, Figure_5 Panel_A-D, Masking_GFF_OUT/produce_masked_fna.py and every RepeatMasker/*.py now take their paths from config.py (no absolute paths left; checked by scan).
Panel_B read a file named ...PC0UW000000.csv; it now reads 4k_kmers_normalised_MCRRC_PC0UW.csv like the other panels.
Panels A-D find the pooled/class CSVs in Figure_5/outputs/ (written by the counting scripts) or, if absent, in Figure_5/ (the committed test file) and write their images to Figure_5/outputs/.
RepeatMasker: Genome_file_dividing reads .fna or .fna.gz, keeps the original unless --delete-original; remove_some_results is a dry run unless --delete; Automation.py fixed (config paths, resume, no crash on the cleanup branch) and stops with a clear message if eukaryotes_for_chargaff.py is missing.
Download scripts write Fasta_GFF_inputs/<Group>/ncbi_dataset/data/... (Viridiplantae -> Plantae); optional RefSeq FTP files go to Fasta_GFF_inputs/RefSeq_FTP/.
New: Download_Data/README.md (NCBI Datasets on Ubuntu), RepeatMasker/README.md, docs/PIPELINE.md; updated README.md, Dataframes/README.md, Fasta_GFF_inputs/README.md, .gitignore.

Tested: all 63 scripts compile; rescan shows 0 absolute paths and 0 undefined config names; Figure 3 runs on your 120 Mammals_CSV folders; Automation/split/transfer/progress/cleanup scripts and the download scripts' syntax tested on synthetic data with a stub RepeatMasker.
NOT tested: real RepeatMasker, real NCBI downloads, produce_masked_fna.py and Figure_5 Panel A-D on the real pooled CSV (not in the zip).
