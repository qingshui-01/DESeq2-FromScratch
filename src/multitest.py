"""M4：Benjamini-Hochberg 多重检验校正（FDR 控制）。

要解决的问题
------------
同时检验两万多个基因时，即使全部是噪声，按 p < 0.05 也会有约 1000 个「显著」结果。
BH 不追求「一个假阳性都没有」，而是把「报出来的显著结果里假阳性所占比例」控制在 5% 以内。

算法（与 R 的 p.adjust(p, method="BH") 完全一致，也是 DESeq2 内部使用的实现）
--------------------------------------------------------------------------
设非缺失 p 值共 n 个，按升序排列为 p(1) <= p(2) <= ... <= p(n)，秩为 i：

    1. 对每个秩计算       q(i) = p(i) * n / i
    2. 从最大的 i 往小扫，保证单调不降：  q(i) = min(q(i), q(i+1))
       （等价于 R 的 cummin(n/i * p[o])）
    3. 截断在 1 以内：    q(i) = min(q(i), 1)
    4. 缺失的 p 值 -> 缺失的 padj

第 2 步（单调化）最容易漏。漏了就会出现「p 值更小的基因反而 padj 更大」这种
数学上不可能的结果，所以单元测试里专门盯这一条。

n 的取值
-------
与 R 一致：n = 传入的非缺失 p 值个数。所以当上游做了独立过滤（independent
filtering）时，只有通过过滤的基因参与校正，这正是本项目要额外导出
deseq2_result_noFilter.csv 的原因（见 REPORT.md 的说明）。
"""


def bh_adjust(pvalues):
    """对一组 p 值做 Benjamini-Hochberg 校正。

    参数
    ----
    pvalues : 序列，元素为 float 或 None（None/NaN 表示缺失）

    返回
    ----
    与输入等长的列表；有效位置为校正后的 padj（float），缺失位置为 None。
    """
    n_total = len(pvalues)

    # 1) 收集有效的 (原始下标, p 值)
    valid = []
    for idx, p in enumerate(pvalues):
        if p is None:
            continue
        p = float(p)
        if p != p:  # NaN 判定（NaN != NaN）
            continue
        # 与 R 一致：p 值先夹到 [0, 1]
        if p < 0.0:
            p = 0.0
        elif p > 1.0:
            p = 1.0
        valid.append((idx, p))

    result = [None] * n_total
    n = len(valid)
    if n == 0:
        return result

    # 2) 按 p 值升序排序（稳定排序，保证相同 p 值顺序可复现）
    valid.sort(key=lambda pair: pair[1])

    # 3) 计算 n/i * p(i)，并从大到小做 cummin
    adjusted = [0.0] * n
    running_min = float("inf")
    for i in range(n - 1, -1, -1):          # i 为 0 基下标，秩 = i + 1
        rank = i + 1
        q = valid[i][1] * n / rank
        if q > 1.0:
            q = 1.0
        if q < running_min:
            running_min = q
        adjusted[i] = running_min

    # 4) 写回到原来的位置
    for i in range(n):
        result[valid[i][0]] = adjusted[i]

    return result
