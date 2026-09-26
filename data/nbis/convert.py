"""把 NBISweden/workshop-RNAseq 的 gene_counts.csv 转成
本项目的 counts.csv / design.csv 格式（基因 ID 列 + 各样本计数列；
design 两列：sample_id, condition，day00=0 / day07=1）。"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "gene_counts.csv")

with open(SRC, newline="", encoding="utf-8-sig") as f:
    r = csv.reader(f)
    header = next(r)
    samples = [h.strip().strip('"') for h in header[1:]]
    rows = []
    for row in r:
        if not row:
            continue
        g = row[0].strip().strip('"')
        vals = [int(round(float(v))) for v in row[1:]]
        rows.append([g] + vals)

with open(os.path.join(HERE, "counts.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["gene_id"] + samples)
    w.writerows(rows)

with open(os.path.join(HERE, "design.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["sample", "group"])
    for s in samples:
        w.writerow([s, 0 if s.startswith("DSSd00") else 1])

print("counts: %d genes x %d samples" % (len(rows), len(samples)))
print("samples:", samples)
