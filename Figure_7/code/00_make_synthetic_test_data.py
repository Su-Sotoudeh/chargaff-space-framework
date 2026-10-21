"""
00_make_synthetic_test_data.py
-------------------------------
Synthetic dataset matching the real schema INCLUDING Organ_ID / Chr_ID, used
only to check that the pipeline runs end-to-end (and that the organism-level
split really keeps organisms apart). NOT real biology.

Structure mimics the real data:
  * every organism has an organism-level GC and strand-composition offset
    shared by ALL its chromosomes (this is what makes a chromosome-level
    split leak);
  * chromosomes of one organism differ by a smaller chromosome-level noise;
  * whole-chromosome composition = Z-weighted mix of the +annotation and the
    -annotation (reverse-complement-like) composition + unannotated noise;
  * a few rows have Z = 0 or Z = 100 (formula methods cannot be evaluated);
  * one virus subgroup has only 6 organisms (ML must skip it, formulas not).
"""
import numpy as np
import pandas as pd
import importlib

rng = np.random.default_rng(0)

# group: (n_organisms, (min_chr, max_chr), organism GC sd, chromosome noise)
SPEC = {
    "Chordata": (40, (8, 16), 0.03, 0.8),
    "Arthropoda": (40, (3, 8), 0.04, 0.9),
    "Plantae": (30, (5, 10), 0.03, 0.8),
    "Fungi": (40, (3, 9), 0.04, 0.9),
    "Archaea": (30, (1, 2), 0.08, 1.0),
    "Bacteria": (80, (1, 3), 0.09, 1.0),
    "Chloroplast": (50, (1, 1), 0.04, 0.8),
    "Mitochondria": (60, (1, 1), 0.07, 1.2),
    "Plasmids": (60, (1, 6), 0.08, 1.5),
}
VIRUS = {("DNA", "+"): 30, ("DNA", "-"): 25, ("DNA", "+/-"): 30,
         ("RNA", "+"): 30, ("RNA", "-"): 6, ("RNA", "+/-"): 25}


def make_rows(group, genome_type, coding, n_org, chr_range, gc_sd, chr_noise, prefix):
    rows = []
    for o in range(n_org):
        oid = f"{prefix}_{o:04d}"
        gc = np.clip(rng.normal(0.45, gc_sd * 2), 0.2, 0.75)
        d = rng.normal(0, 0.02)       # A-T asymmetry of the sense composition
        e = rng.normal(0, 0.02)       # G-C asymmetry
        z_org = np.clip(rng.normal(55, 12), 15, 90)
        for k in range(rng.integers(chr_range[0], chr_range[1] + 1)):
            noise = rng.normal(0, chr_noise / 100, 4)
            a = 0.5 * (1 - gc) + d + noise[0]
            t = 0.5 * (1 - gc) - d + noise[1]
            g = 0.5 * gc + e + noise[2]
            c = 0.5 * gc - e + noise[3]
            s = a + t + g + c
            ap, tp, gp, cp = 100 * a / s, 100 * t / s, 100 * g / s, 100 * c / s
            z = float(np.clip(z_org + rng.normal(0, 8), 3, 97))
            u = rng.random()
            if u < 0.01:
                z = 0.0
            elif u < 0.02:
                z = 100.0
            # -annotation as seen on the reference strand ~ complement of +annotation
            cn = rng.normal(0, 1.0, 4)
            ac, tc, gc_, cc = tp + cn[0], ap + cn[1], cp + cn[2], gp + cn[3]
            s2 = ac + tc + gc_ + cc
            ac, tc, gc_, cc = 100 * ac / s2, 100 * tc / s2, 100 * gc_ / s2, 100 * cc / s2
            un = rng.normal(0, 0.4, 4)          # unannotated sequence noise
            aw = (z * ap + (100 - z) * ac) / 100 + un[0]
            tw = (z * tp + (100 - z) * tc) / 100 + un[1]
            gw = (z * gp + (100 - z) * gc_) / 100 + un[2]
            cw = (z * cp + (100 - z) * cc) / 100 + un[3]
            sw = aw + tw + gw + cw
            aw, tw, gw, cw = 100 * aw / sw, 100 * tw / sw, 100 * gw / sw, 100 * cw / sw
            rows.append(dict(
                Organ_ID=oid, Chr_ID=f"{oid}_chr{k + 1}", Group=group,
                Genome_Type=genome_type, Coding=coding,
                A_W=aw, T_W=tw, G_W=gw, C_W=cw,
                Z=z, P_portion=z, C_portion=100 - z,
                annot_P_A=ap, annot_P_T=tp, annot_P_G=gp, annot_P_C=cp,
                annot_C_A=ac, annot_C_T=tc, annot_C_G=gc_, annot_C_C=cc))
    return rows


rows = []
for g, (n_org, chr_range, gc_sd, chr_noise) in SPEC.items():
    rows += make_rows(g, "DNA", "+/-", n_org, chr_range, gc_sd, chr_noise, g[:4].upper())
for (gt, cd), n_org in VIRUS.items():
    rows += make_rows("Viruses", gt, cd, n_org, (1, 1), 0.10, 2.0, f"V{gt}{cd}".replace("/", "").replace("+", "p").replace("-", "m"))

df = pd.DataFrame(rows)
paths = importlib.import_module("paths")
out_path = paths.table_path("synthetic_test_data.csv")
df.to_csv(out_path, index=False)
print(f"Synthetic test data created: {df.shape}  ->  {out_path}")
print(df.groupby("Group").agg(rows=("Chr_ID", "size"), organisms=("Organ_ID", "nunique")))
