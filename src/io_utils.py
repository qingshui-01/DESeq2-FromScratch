"""输入/输出工具：读取计数矩阵与分组表，写出结果表。

输入格式约定（对拍的两条路径必须吃完全相同的文件）
--------------------------------------------------
counts.csv
    第 1 行是表头：第一个字段为空（或行名占位），后面是样本名
    之后每行一个基因：第一个字段是基因 ID，后面是各样本的计数（非负整数）

    ,SRR1039508,SRR1039509,...
    ENSG00000000003,679,448,...
    ENSG00000000005,0,0,...

design.csv
    表头固定为 sample,condition,cell,group
    group 是 0/1，1 表示处理组（log2FoldChange 的分子）

结果表 result.csv
    表头与 DESeq2 results() 的输出列名完全对齐，便于逐列对拍：
    gene_id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj

编码统一 UTF-8，写出时不加 BOM。
"""

import csv


def _strip_quotes(s):
    if s is None:
        return s
    s = s.strip()
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        return s[1:-1]
    return s


def read_counts_csv(path):
    """读取计数矩阵。

    返回 (gene_ids, sample_ids, counts)，其中 counts[基因][样本] 为 int。
    """
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("counts.csv 是空文件")
        sample_ids = [_strip_quotes(x) for x in header[1:]]
        if not sample_ids:
            raise ValueError("counts.csv 表头里没有样本列")

        gene_ids = []
        counts = []
        n_samples = len(sample_ids)
        for line_no, row in enumerate(reader, start=2):
            if not row or all(_strip_quotes(c) == "" for c in row):
                continue
            if len(row) != n_samples + 1:
                raise ValueError(
                    "counts.csv 第 %d 行有 %d 个字段，期望 %d 个"
                    % (line_no, len(row), n_samples + 1))
            gene_ids.append(_strip_quotes(row[0]))
            vals = []
            for cell in row[1:]:
                cell = _strip_quotes(cell)
                try:
                    v = float(cell)
                except ValueError:
                    raise ValueError(
                        "counts.csv 第 %d 行的值不是数字：%r" % (line_no, cell))
                if v < 0:
                    raise ValueError("counts.csv 第 %d 行出现负值" % line_no)
                if v != int(v):
                    raise ValueError(
                        "counts.csv 第 %d 行出现非整数计数（本项目要求原始计数）"
                        % line_no)
                vals.append(int(v))
            counts.append(vals)

    if not counts:
        raise ValueError("counts.csv 没有数据行")
    return gene_ids, sample_ids, counts


def read_design_csv(path):
    """读取分组表。返回 dict(样本名 -> dict(condition, cell, group))。"""
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError("design.csv 是空文件")
        fields = [_strip_quotes(x) for x in reader.fieldnames]
        if "sample" not in fields or "group" not in fields:
            raise ValueError(
                "design.csv 必须包含 sample 与 group 两列，实际为：%s" % fields)
        out = {}
        for row in reader:
            rec = {_strip_quotes(k): _strip_quotes(v)
                   for k, v in row.items() if k is not None}
            sample = rec.get("sample")
            if not sample:
                continue
            try:
                group = int(float(rec["group"]))
            except (TypeError, ValueError):
                raise ValueError("design.csv 中样本 %s 的 group 不是数字" % sample)
            if group not in (0, 1):
                raise ValueError(
                    "design.csv 中样本 %s 的 group = %d，本项目只支持 0/1 两组"
                    % (sample, group))
            out[sample] = {
                "condition": rec.get("condition", ""),
                "cell": rec.get("cell", ""),
                "group": group,
            }
    if not out:
        raise ValueError("design.csv 没有数据行")
    return out


def align_design(sample_ids, design):
    """把 counts 的样本顺序与 design 对齐，返回 (groups, conditions, cells)。

    缺失样本或分组不全是 0/1 都会明确报错，因为对拍的前提是两边样本完全一致。
    """
    groups, conditions, cells = [], [], []
    for s in sample_ids:
        if s not in design:
            raise ValueError("样本 %s 出现在 counts.csv 但不在 design.csv 中" % s)
        groups.append(design[s]["group"])
        conditions.append(design[s]["condition"])
        cells.append(design[s]["cell"])
    if len(set(groups)) != 2:
        raise ValueError(
            "design 中只有一组（group 取值 %s），本项目要求恰好两组"
            % sorted(set(groups)))
    return groups, conditions, cells


def write_results_csv(path, gene_ids, rows):
    """写出结果表。

    rows 是与 gene_ids 等长的列表，每项是一个 dict，键为
    baseMean / log2FoldChange / lfcSE / stat / pvalue / padj（值可为 None）。
    """
    header = ["gene_id", "baseMean", "log2FoldChange", "lfcSE",
              "stat", "pvalue", "padj"]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for gid, row in zip(gene_ids, rows):
            writer.writerow([gid] + [_fmt(row.get(k)) for k in header[1:]])
    return path


def _fmt(v):
    """数值格式化：None/NaN 写成空字符串，其余保留足够有效位。"""
    if v is None:
        return ""
    if v != v:      # NaN
        return ""
    if isinstance(v, bool):
        return "TRUE" if v else "FALSE"
    if v == int(v) and abs(v) < 1e15:
        return str(int(v))
    return repr(float(v))


def read_results_csv(path):
    """读回结果表（用于对拍）。返回 dict(gene_id -> dict(列名 -> float/None))。"""
    out = {}
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rec = {}
            gid = None
            for k, v in row.items():
                if k is None:
                    continue
                key = _strip_quotes(k)
                val = _strip_quotes(v) if v is not None else ""
                if key == "gene_id":
                    gid = val
                    continue
                if val == "" or val.upper() in ("NA", "NAN", "NULL"):
                    rec[key] = None
                else:
                    try:
                        rec[key] = float(val)
                    except ValueError:
                        rec[key] = None
            if gid is not None:
                out[gid] = rec
    return out


def write_generic_csv(path, header, rows):
    """写任意 CSV（用于导出中间量）。"""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for row in rows:
            writer.writerow([_fmt(v) if isinstance(v, float) else
                             ("" if v is None else v) for v in row])
    return path
