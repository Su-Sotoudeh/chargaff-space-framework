import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(next(p for p in _Path(__file__).resolve().parents if (p / "config.py").exists())))
import config  # noqa: E402  (all input/output paths live in config.py)
import os
import re
import gzip
import time
import collections
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from Bio import SeqIO


# ----------------- Configuration -----------------
MIN_CHROMOSOME_BP = 5_000_000
MAX_WORKERS = 8

INPUT_PATH = str(config.require(config.MASKED_FASTAS, "Put the masked FASTA files into Fasta_GFF_inputs/Masked_Fastas/<organism>/ (see Fasta_GFF_inputs/README.md)."))   # Fasta_GFF_inputs/Masked_Fastas/<organism>/
GROUP = "Mammals"


# ----------------- Symbol Sets -----------------
CDS_P_SYM = set("ZXVB")
CDS_C_SYM = set("zxvb")

TE_P_SYM = set(["1", "2", "3", "4", "L", "M", "@", "#", "O", "P", "J", "K", "*", "Y", "U", "I"])
TE_C_SYM = set(["5", "6", "7", "8", "l", "m", "$", "%", "o", "p", "j", "k", "+", "y", "u", "i"])


# ----------------- Helpers -----------------
def safe_div(a, b, default=0):
    if b == 0:
        return default
    return round(a / b, 4)


def is_valid_chromosome_name(description, seq_len):
    if seq_len < MIN_CHROMOSOME_BP:
        return False

    text = re.sub(r"\s+", " ", str(description).lower()).strip()
    bad_keywords = [
        "scaffold",
        "contig",
        "unlocalized",
        "unplaced",
        "patch",
        "plasmid",
        "mitochond",
        "chloroplast",
    ]
    if any(word in text for word in bad_keywords):
        return False

    if re.search(r"\bchr[\s._-]*[0-9xyzwm]+\b", text):
        return True
    if str(description).upper().startswith("NC_"):
        return True
    if "chromosome" in text or "complete sequence" in text or "complete genome" in text:
        return True

    return False


def count_segments(sequence, symbols_set):
    in_segment = False
    count = 0

    for ch in sequence:
        if ch in symbols_set:
            if not in_segment:
                count += 1
                in_segment = True
        else:
            in_segment = False

    return count


def nucleotide_analysis_decoded(sequence):
    """
    Rows:
    +Int, -Int, +UTR, -UTR, +CDS, -CDS, +TEs, -TEs, Shared, Unannotated
    Columns:
    A, T, C, G
    """
    data = np.zeros((10, 4), dtype=np.int64)

    sym_map = {
        # Int +
        "Q": (0, 0), "W": (0, 1), "E": (0, 2), "R": (0, 3),
        # Int -
        "q": (1, 0), "w": (1, 1), "e": (1, 2), "r": (1, 3),

        # UTR +
        "S": (2, 0), "D": (2, 1), "F": (2, 2), "H": (2, 3),
        # UTR -
        "s": (3, 0), "d": (3, 1), "f": (3, 2), "h": (3, 3),

        # CDS +
        "Z": (4, 0), "X": (4, 1), "V": (4, 2), "B": (4, 3),
        # CDS -
        "z": (5, 0), "x": (5, 1), "v": (5, 2), "b": (5, 3),

        # TEs +
        "1": (6, 0), "2": (6, 1), "3": (6, 2), "4": (6, 3),
        "L": (6, 0), "M": (6, 1), "@": (6, 2), "#": (6, 3),
        "O": (6, 0), "P": (6, 1), "J": (6, 2), "K": (6, 3),
        "*": (6, 0), "Y": (6, 1), "U": (6, 2), "I": (6, 3),

        # TEs -
        "5": (7, 0), "6": (7, 1), "7": (7, 2), "8": (7, 3),
        "l": (7, 0), "m": (7, 1), "$": (7, 2), "%": (7, 3),
        "o": (7, 0), "p": (7, 1), "j": (7, 2), "k": (7, 3),
        "+": (7, 0), "y": (7, 1), "u": (7, 2), "i": (7, 3),

        # Shared
        "A": (8, 0), "T": (8, 1), "C": (8, 2), "G": (8, 3),

        # Unannotated
        "a": (9, 0), "t": (9, 1), "c": (9, 2), "g": (9, 3),
    }

    counts = collections.Counter(sequence)
    for sym, count in counts.items():
        if sym in sym_map:
            row, col = sym_map[sym]
            data[row, col] += count

    df = pd.DataFrame(data, columns=["A", "T", "C", "G"])
    df.index = [
        "+Int", "-Int",
        "+UTR", "-UTR",
        "+CDS", "-CDS",
        "+TEs", "-TEs",
        "Shared", "Unannotated",
    ]
    return df


def calc_element_stats(a, t, c, g, length):
    total = a + t + c + g
    if total == 0:
        return {
            "GC": 0,
            "portion": 0,
            "A_ano": 0,
            "T_ano": 0,
            "C_ano": 0,
            "G_ano": 0,
            "A2T": 0,
            "G2C": 0,
            "count": 0,
        }

    return {
        "GC": round(((c + g) / total) * 100, 4),
        "portion": round((total / length) * 100, 4),
        "A_ano": round((a / total) * 100, 4),
        "T_ano": round((t / total) * 100, 4),
        "C_ano": round((c / total) * 100, 4),
        "G_ano": round((g / total) * 100, 4),
        "A2T": safe_div(a, t),
        "G2C": safe_div(g, c),
        "count": total,
    }


def analyze_one_element(element_name, plus_counts, minus_counts, length):
    ap, tp, cp, gp = plus_counts
    ac, tc, cc, gc = minus_counts

    p_stats = calc_element_stats(ap, tp, cp, gp, length)
    c_stats = calc_element_stats(ac, tc, cc, gc, length)

    smaller = min(p_stats["portion"], c_stats["portion"])
    bigger = max(p_stats["portion"], c_stats["portion"])
    s2b = round(smaller / bigger, 4) if bigger != 0 else 0

    total_pc = p_stats["count"] + c_stats["count"]
    dif2sum = round(abs(p_stats["count"] - c_stats["count"]) / (2 * total_pc) * 100, 4) if total_pc != 0 else 0

    return {
        f"{element_name}_GC_P": p_stats["GC"],
        f"{element_name}_GC_C": c_stats["GC"],
        f"{element_name}_P%": p_stats["portion"],
        f"{element_name}_C%": c_stats["portion"],
        f"{element_name}_P2C": safe_div(p_stats["count"], c_stats["count"]),
        f"{element_name}_S2B": s2b,
        f"{element_name}_dif2sum": dif2sum,
        f"{element_name}_A2T_ano": safe_div(ap + ac, tp + tc),
        f"{element_name}_G2C_ano": safe_div(gp + gc, cp + cc),
        f"{element_name}_P_A_ano": p_stats["A_ano"],
        f"{element_name}_P_T_ano": p_stats["T_ano"],
        f"{element_name}_P_C_ano": p_stats["C_ano"],
        f"{element_name}_P_G_ano": p_stats["G_ano"],
        f"{element_name}_P_A2T": p_stats["A2T"],
        f"{element_name}_P_G2C": p_stats["G2C"],
        f"{element_name}_C_A_ano": c_stats["A_ano"],
        f"{element_name}_C_T_ano": c_stats["T_ano"],
        f"{element_name}_C_C_ano": c_stats["C_ano"],
        f"{element_name}_C_G_ano": c_stats["G_ano"],
        f"{element_name}_C_A2T": c_stats["A2T"],
        f"{element_name}_C_G2C": c_stats["G2C"],
    }


def build_column_order():
    cols = [
        "Name",
        "Chr_ID",
        "Length(M)",
        "GC_W",
        "Wh_ano",
        "Int_ano",
        "UTR_ano",
        "CDS_ano",
        "TEs_ano",
        "share_ano",
        "Un_ano",
        "sum",
        "A2T_W",
        "G2C_W",
        "TEs_P_count",
        "TEs_C_count",
        "CDS_P_count",
        "CDS_C_count",
    ]

    elements = ["Int", "UTR", "CDS", "TEs", "Share", "Unanot"]
    for el in elements:
        cols.extend([
            f"{el}_GC_P",
            f"{el}_GC_C",
            f"{el}_P%",
            f"{el}_C%",
            f"{el}_P2C",
            f"{el}_S2B",
            f"{el}_dif2sum",
            f"{el}_A2T_ano",
            f"{el}_G2C_ano",
            f"{el}_P_A_ano",
            f"{el}_P_T_ano",
            f"{el}_P_C_ano",
            f"{el}_P_G_ano",
            f"{el}_P_A2T",
            f"{el}_P_G2C",
            f"{el}_C_A_ano",
            f"{el}_C_T_ano",
            f"{el}_C_C_ano",
            f"{el}_C_G_ano",
            f"{el}_C_A2T",
            f"{el}_C_G2C",
        ])

    return cols


COLUMN_ORDER = build_column_order()


# ----------------- Worker -----------------
def process_one_organism(organism, input_root):
    try:
        start_time = time.time()
        print(f"[START] {organism}")

        org_dir = os.path.join(input_root, organism)
        fasta_files = [
            f for f in os.listdir(org_dir)
            if f.lower().endswith((".fna.gz", ".fa.gz", ".fasta.gz"))
        ]

        if not fasta_files:
            print(f"[ERROR] {organism}: no fasta file found")
            return None

        organism_rows = []

        for fasta_name in fasta_files:
            fasta_path = os.path.join(org_dir, fasta_name)

            with gzip.open(fasta_path, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    seq = str(record.seq)

                    if not is_valid_chromosome_name(record.description, len(seq)):
                        continue

                    length = len(seq) - seq.count("n") - seq.count("N")
                    if length <= 0:
                        continue

                    te_p_count = count_segments(seq, TE_P_SYM)
                    te_c_count = count_segments(seq, TE_C_SYM)
                    cds_p_count = count_segments(seq, CDS_P_SYM)
                    cds_c_count = count_segments(seq, CDS_C_SYM)

                    df = nucleotide_analysis_decoded(seq)

                    a_w = int(df["A"].sum())
                    t_w = int(df["T"].sum())
                    c_w = int(df["C"].sum())
                    g_w = int(df["G"].sum())
                    total_w = a_w + t_w + c_w + g_w

                    if total_w == 0:
                        continue

                    main_dict = {
                        "Name": organism,
                        "Chr_ID": record.id,
                        "Length(M)": round(length / 1_000_000, 6),
                        "GC_W": round(((g_w + c_w) / total_w) * 100, 4),
                        "A2T_W": safe_div(a_w, t_w),
                        "G2C_W": safe_div(g_w, c_w),
                        "TEs_P_count": te_p_count,
                        "TEs_C_count": te_c_count,
                        "CDS_P_count": cds_p_count,
                        "CDS_C_count": cds_c_count,
                    }

                    int_plus = tuple(df.loc["+Int", ["A", "T", "C", "G"]].tolist())
                    int_minus = tuple(df.loc["-Int", ["A", "T", "C", "G"]].tolist())
                    main_dict.update(analyze_one_element("Int", int_plus, int_minus, length))

                    utr_plus = tuple(df.loc["+UTR", ["A", "T", "C", "G"]].tolist())
                    utr_minus = tuple(df.loc["-UTR", ["A", "T", "C", "G"]].tolist())
                    main_dict.update(analyze_one_element("UTR", utr_plus, utr_minus, length))

                    cds_plus = tuple(df.loc["+CDS", ["A", "T", "C", "G"]].tolist())
                    cds_minus = tuple(df.loc["-CDS", ["A", "T", "C", "G"]].tolist())
                    main_dict.update(analyze_one_element("CDS", cds_plus, cds_minus, length))

                    tes_plus = tuple(df.loc["+TEs", ["A", "T", "C", "G"]].tolist())
                    tes_minus = tuple(df.loc["-TEs", ["A", "T", "C", "G"]].tolist())
                    main_dict.update(analyze_one_element("TEs", tes_plus, tes_minus, length))

                    shared_row = df.loc["Shared", ["A", "T", "C", "G"]].tolist()
                    shared_plus = tuple(x / 2 for x in shared_row)
                    shared_minus = tuple(x / 2 for x in shared_row)
                    main_dict.update(analyze_one_element("Share", shared_plus, shared_minus, length))

                    unanot_row = df.loc["Unannotated", ["A", "T", "C", "G"]].tolist()
                    unanot_plus = tuple(x / 2 for x in unanot_row)
                    unanot_minus = tuple(x / 2 for x in unanot_row)
                    main_dict.update(analyze_one_element("Unanot", unanot_plus, unanot_minus, length))

                    int_ano = round((sum(int_plus) + sum(int_minus)) / length * 100, 4)
                    utr_ano = round((sum(utr_plus) + sum(utr_minus)) / length * 100, 4)
                    cds_ano = round((sum(cds_plus) + sum(cds_minus)) / length * 100, 4)
                    tes_ano = round((sum(tes_plus) + sum(tes_minus)) / length * 100, 4)
                    share_ano = round((sum(shared_plus) + sum(shared_minus)) / length * 100, 4)
                    un_ano = round((sum(unanot_plus) + sum(unanot_minus)) / length * 100, 4)

                    wh_ano = round(int_ano + utr_ano + cds_ano + tes_ano, 4)
                    total_sum = round(wh_ano + share_ano + un_ano, 4)

                    res_df = pd.DataFrame([main_dict])

                    res_df.insert(4, "Wh_ano", wh_ano)
                    res_df.insert(5, "Int_ano", int_ano)
                    res_df.insert(6, "UTR_ano", utr_ano)
                    res_df.insert(7, "CDS_ano", cds_ano)
                    res_df.insert(8, "TEs_ano", tes_ano)
                    res_df.insert(9, "share_ano", share_ano)
                    res_df.insert(10, "Un_ano", un_ano)
                    res_df.insert(11, "sum", total_sum)

                    for col in COLUMN_ORDER:
                        if col not in res_df.columns:
                            res_df[col] = 0

                    res_df = res_df.reindex(columns=COLUMN_ORDER)
                    organism_rows.append(res_df)

        if organism_rows:
            out_df = pd.concat(organism_rows, ignore_index=True)
            elapsed = time.time() - start_time
            print(f"[DONE] {organism} | rows={len(out_df)} | time={elapsed:.2f}s")
            return out_df

        print(f"[DONE] {organism} | no valid chromosomes found")
        return None

    except Exception as e:
        print(f"[ERROR] {organism}: {e}")
        return None


# ----------------- Main Runner -----------------
def run_analysis(input_root, group_name):
    out_dir = os.path.join(str(config.output_dir("GFF_Mapping", group_name)), f"{group_name}_Analysis_{datetime.now().strftime('%Y-%m-%d')}")
    os.makedirs(out_dir, exist_ok=True)

    master_csv = os.path.join(out_dir, f"{group_name}_fixed.csv")

    all_orgs = sorted([
        d for d in os.listdir(input_root)
        if os.path.isdir(os.path.join(input_root, d))
    ])

    if os.path.exists(master_csv):
        try:
            processed = set(pd.read_csv(master_csv, usecols=["Name"])["Name"].dropna().unique())
        except Exception:
            processed = set()
    else:
        processed = set()

    to_process = [org for org in all_orgs if org not in processed]

    print(f"[INIT] Total: {len(all_orgs)} | Processed: {len(processed)} | To Process: {len(to_process)}")

    if not to_process:
        print("[PROGRESS] Nothing left to process")
        return

    completed = 0

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_map = {
            executor.submit(process_one_organism, org, input_root): org
            for org in to_process
        }

        for future in as_completed(future_map):
            organism = future_map[future]
            result = future.result()

            if result is not None and not result.empty:
                write_header = not os.path.exists(master_csv)
                result.to_csv(master_csv, mode="a", index=False, header=write_header)
                print(f"[SAVED] {organism}")

            completed += 1
            print(f"[PROGRESS] {completed}/{len(to_process)} finished")


if __name__ == "__main__":
    run_analysis(INPUT_PATH, GROUP)
