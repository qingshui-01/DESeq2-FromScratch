# -*- coding: utf-8 -*-
"""反向验证（故障注入）：把已知 bug 重新注入，确认验收点会变红。

这么做的理由
------------
「验收点全绿」本身不能说明验收点有用：**一条永远为真的断言也是绿的**。
唯一能证明验收点有用的办法是：把已知 bug 放回去，看它会不会变红。

本脚本对每个故障都是「备份 → 注入 → 跑验收 → 还原（finally 保证）」，
原件不会被留下改动。

用法
----
    python tests/fault_injection.py            # 只跑代码级注入（约 30 秒）
    python tests/fault_injection.py --e2e      # 额外做一次端到端注入：
                                               # 注入后**重跑 deg.py**，看
                                               # I-1/I-2 会不会变红（约 4 分钟）

--e2e 为什么必要：I-1/I-2/I-4 读的是**已产出的** data/my_result.csv。
只改 src/ 而不重跑 deg.py，它们仍然全绿，那是在检查缓存，不是在检查代码。
"""

import argparse
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BACKUP = os.path.join(ROOT, "build", "fault_injection_backup")

# name -> (相对路径, 原文, 注入后, 期望变红的验收点, 是否需要重跑 deg.py)
FAULTS = {
    "1 工作响应未用截断后的 mu": (
        "src/glm_nb.py",
        "        z = [(math.log(mu[j]) - log_sf[j]) + (counts[j] - mu[j]) / mu[j]\n"
        "             for j in range(m)]",
        "        z = [eta[j] + log_sf[j] + (counts[j] - mu[j]) / mu[j]\n"
        "             for j in range(m)]",
        "F-3 / F-5 / I-1 / I-2",
        False,
    ),
    "2 收敛判据比较顺序写反": (
        "src/glm_nb.py",
        "        delta_beta = max(abs(beta_try[a] - beta[a]) for a in range(p))\n"
        "        beta = beta_try",
        "        beta = beta_try\n"
        "        delta_beta = max(abs(beta_try[a] - beta[a]) for a in range(p))",
        "F-5",
        False,
    ),
    "3 离散度 gene-wise 缺上界": (
        "src/dispersion.py",
        "        v = max(v, MIN_DISP)\n        v = min(v, max_disp)",
        "        v = max(v, MIN_DISP)",
        "E-5",
        False,
    ),
    "4 丢失 Cox-Reid 校正项": (
        "src/dispersion.py",
        "        v = gene_wise_dispersion(counts[i], mu, alpha_start[i],\n"
        "                                 X=X, use_cr=True, initial_lp=initial_ll)",
        "        v = gene_wise_dispersion(counts[i], mu, alpha_start[i],\n"
        "                                 X=X, use_cr=False, initial_lp=initial_ll)",
        "I-5",
        True,   # I-5 读 dispersion_mine.csv，必须重跑 deg.py 才会变红
    ),
}


def paths(rel):
    p = os.path.join(ROOT, rel.replace("/", os.sep))
    return p, os.path.join(BACKUP, rel.replace("/", "_"))


def backup_all():
    os.makedirs(BACKUP, exist_ok=True)
    for rel in {f[0] for f in FAULTS.values()}:
        p, b = paths(rel)
        shutil.copyfile(p, b)


def inject(rel, old, new):
    p, _b = paths(rel)
    with open(p, "r", encoding="utf-8") as fh:
        src = fh.read()
    if src.count(old) != 1:
        raise SystemExit("注入失败：原文出现 %d 次（应为 1）\n%s"
                         % (src.count(old), old))
    with open(p, "w", encoding="utf-8", newline="") as fh:
        fh.write(src.replace(old, new))


def restore_all():
    for rel in {f[0] for f in FAULTS.values()}:
        p, b = paths(rel)
        shutil.copyfile(b, p)


def run(cmd, timeout=1200):
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                          encoding="utf-8", errors="replace", timeout=timeout)


def acceptance_failures():
    p = run([sys.executable, "tests/run_acceptance.py"])
    lines = p.stdout.splitlines()
    failed = [ln.strip() for ln in lines if "[失败]" in ln]
    detail = {}
    for i, ln in enumerate(lines):
        if "[失败]" in ln and i + 1 < len(lines):
            detail[ln.split()[1]] = lines[i + 1].strip()
    summary = next((ln for ln in lines if "汇总：" in ln), "")
    return failed, detail, summary


def main():
    ap = argparse.ArgumentParser(description="故障注入反向验证")
    ap.add_argument("--e2e", action="store_true",
                    help="额外做端到端注入（重跑 deg.py）")
    args = ap.parse_args()

    backup_all()
    ok = True
    try:
        print("=" * 76)
        print("反向验证：注入已知 bug，确认验收点会变红")
        print("=" * 76)
        for name, (rel, old, new, expect, needs_e2e) in FAULTS.items():
            print("\n>>> 注入：%s   (%s)" % (name, rel))
            print("    期望变红：%s" % expect)
            if needs_e2e:
                print("    （该故障需重跑 deg.py 才能被验收点看到，"
                      "仅在 --e2e 模式下注入）")
                continue
            try:
                inject(rel, old, new)
                failed, detail, summary = acceptance_failures()
                print("    %s" % summary)
                for ln in failed:
                    code = ln.split()[1]
                    print("      %s   %s" % (ln, detail.get(code, "")))
                if failed:
                    print("    >>> 验收变红 ✔（这个 bug 被守住了）")
                else:
                    print("    >>> 验收仍然全绿 ✘（存在验收盲区）")
                    ok = False
            finally:
                restore_all()
                print("    （已还原 %s）" % rel)

        if args.e2e:
            # 端到端注入：① 故障 1 证明 I-1/I-2 检查的是代码而不是缓存；
            # ② 所有 needs_e2e 的故障（如丢失 CR 校正 → I-5 应变红）。
            e2e_faults = [list(FAULTS.items())[0]]
            e2e_faults += [(n, f) for n, f in FAULTS.items() if f[4]
                           and n != e2e_faults[0][0]]
            for name, (rel, old, new, expect, _e) in e2e_faults:
                print("\n" + "=" * 76)
                print("端到端注入：%s → **重跑 deg.py** → 看 %s 是否变红"
                      % (name, expect))
                print("=" * 76)
                try:
                    print("\n[0] 注入前（应为绿）")
                    failed, _d, summary = acceptance_failures()
                    print("    %s" % summary)
                    inject(rel, old, new)
                    print("\n[1] 注入后重跑 deg.py（约 90 秒）...")
                    r = run([sys.executable, "deg.py", "--counts",
                             "data/counts.csv", "--design", "data/design.csv",
                             "--out", "data/my_result.csv",
                             "--dump-intermediates", "data/intermediates",
                             "--quiet"])
                    if r.returncode != 0:
                        print("    deg.py 失败：%s" % r.stderr[-800:])
                    failed, detail, summary = acceptance_failures()
                    print("    %s" % summary)
                    for ln in failed:
                        print("      %s   %s" % (ln, detail.get(ln.split()[1], "")))
                    if failed:
                        print("    >>> 变红 ✔（说明验收点检查的是代码，不是缓存）")
                    else:
                        print("    >>> 仍然全绿 ✘")
                        ok = False
                finally:
                    restore_all()
                    print("\n[2] 还原源码并重跑 deg.py（恢复正确产物）...")
                    run([sys.executable, "deg.py", "--counts", "data/counts.csv",
                         "--design", "data/design.csv", "--out",
                         "data/my_result.csv", "--dump-intermediates",
                         "data/intermediates", "--quiet"])
                    failed, _d, summary = acceptance_failures()
                    print("    %s" % summary)
                    if failed:
                        print("    >>> 还原后仍失败 ✘，请检查是否真的还原了")
                        ok = False
                    else:
                        print("    >>> 还原后恢复全绿 ✔")
    finally:
        restore_all()

    print("\n" + "=" * 76)
    print("结论：%s" % ("全部故障都被验收点抓住 ✔" if ok else "存在盲区，需要补验收点 ✘"))
    print("=" * 76)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
