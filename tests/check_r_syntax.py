# -*- coding: utf-8 -*-
"""对 R 脚本做一次轻量的语法体检。

用途
----
`laptop/deseq2_baseline.R` 只能在旧笔记本上运行，而主力机没有 R。
这意味着**写错一个括号就要占用用户另一台机器上的一次试错**。
本模块在主力机上尽量提前发现这类低级错误：

  * 括号 `()` `{}` `[]` 是否平衡
  * 引号是否成对闭合
  * 换行符是否是纯 LF（CRLF 粘进 Linux 终端会多出 \\r，导致命令报错）
  * 已登记函数的具名实参是否都是该函数的合法形参
    （这一条是为 `gc(quiet = TRUE)` 这类错误加的：R 的 gc() 没有 quiet 参数，
     写错参数名会直接中止脚本，而括号检查根本看不出来）

它**不是** R 解析器，只能粗筛。真正的验证必须在笔记本上跑一次。

用法
----
    python tests/check_r_syntax.py <脚本路径>

也可被 tests/run_acceptance.py 作为验收点 B-3 调用。
"""

import os
import re
import sys

BACKSLASH = chr(92)

# ---------------------------------------------------------------------------
# 函数调用「参数名白名单」
# ---------------------------------------------------------------------------
# 只登记 deseq2_baseline.R 里真正用到的、且形参名我逐个核对过官方文档的函数。
# 之所以要这张表：括号平衡检查抓不到 gc(quiet = TRUE) 这种错：R 的 gc() 根本
# 没有 quiet 参数（只有 verbose / reset / full），写错参数名 R 会直接抛
# 「参数没有用」并停止执行，而这类错误每次都只能靠用户跑到笔记本上才发现。
#
# WARNING: 往 R 脚本里新增一个函数调用时，必须同步把它加进这张表，
#    否则新调用的参数名不会被检查（漏检，而不是误报）。
SIGNATURES = {
    "gc": {"verbose", "reset", "full"},
    "DESeq": {"object", "test", "fitType", "sfType", "betaPrior",
              "minReplicatesForReplace", "modelMatrixType", "useT", "minmu",
              "parallel", "BPPARAM", "quiet"},
    "DESeqDataSet": {"se", "design", "ignoreRank", "countData", "colData"},
    "results": {"object", "contrast", "name", "lfcThreshold", "altHypothesis",
                "listValues", "cooksCutoff", "independentFiltering", "alpha",
                "pAdjustMethod", "filter", "theta", "format", "saveCols",
                "test", "addMLE", "tidy", "parallel", "BPPARAM", "minmu"},
    "sizeFactors": {"object", "type", "locfunc", "geoMeans", "controlGenes",
                    "normMatrix"},
    "counts": {"object", "normalized", "replaced"},
    "assay": {"x", "i", "withDimnames"},
    "mcols": {"x", "value", "use.names"},
    "readLines": {"con", "n", "ok", "warn", "encoding"},
    "grep": {"pattern", "x", "ignore.case", "perl", "value", "fixed",
             "useBytes", "invert"},
    "write.csv": {"x", "file", "append", "quote", "sep", "eol", "na", "dec",
                  "row.names", "col.names", "qmethod", "fileEncoding"},
    "dir.create": {"path", "showWarnings", "recursive", "mode"},
    "list.files": {"path", "pattern", "all.files", "full.names", "recursive",
                   "ignore.case", "include.dirs", "no.."},
    "factor": {"x", "levels", "labels", "exclude", "ordered", "nmax"},
    "difftime": {"time1", "time2", "tz", "units"},
    "median": {"x", "na.rm"},
    "min": {"na.rm"},
    "max": {"na.rm"},
    "cor": {"x", "y", "use", "method"},
    "round": {"x", "digits"},
    "merge": {"x", "y", "by", "by.x", "by.y", "all", "all.x", "all.y", "sort",
              "suffixes", "no.dups", "incomparables"},
    "rm": {"list", "envir"},
    "library": {"package", "help", "pos", "lib.loc", "character.only",
                "logical.return", "warn.conflicts", "verbose", "quietly"},
}


def _scan_calls(src):
    """扫出源码里所有「已登记函数」的具名实参。

    返回 [(行号, 函数名, [参数名, ...]), ...]。
    只做词法扫描（跳过字符串与 # 注释），并维护一个括号栈来保证
    `f(a = 1, g(b = 2))` 里的 b 不会被算到 f 头上。
    """
    calls = []
    n = len(src)
    i, line = 0, 1
    in_str = esc = in_comment = False
    ident = ""
    stack = []

    while i < n:
        c = src[i]
        if c == "\n":
            line += 1
            in_comment = False
            i += 1
            continue
        if in_comment:
            i += 1
            continue
        if in_str:
            if esc:
                esc = False
            elif c == BACKSLASH:
                esc = True
            elif c == '"':
                in_str = False
            i += 1
            continue
        if c == "#":
            in_comment = True
            i += 1
            continue
        if c == '"':
            in_str = True
            ident = ""
            i += 1
            continue
        if c.isalnum() or c in "._":
            ident += c
            i += 1
            continue
        if c == "(":
            if stack:
                stack[-1]["at"] += 1
            stack.append({"name": ident if ident in SIGNATURES else None,
                          "line": line, "args": [], "at": 0})
            ident = ""
            i += 1
            continue
        if c == ")":
            if stack:
                fr = stack.pop()
                if fr["name"]:
                    calls.append((fr["line"], fr["name"], fr["args"]))
                if stack:
                    stack[-1]["at"] -= 1
            ident = ""
            i += 1
            continue
        if c == "=":
            prev = src[i - 1] if i > 0 else ""
            nxt = src[i + 1] if i + 1 < n else ""
            # 排除 ==、<=、>=、!=
            if nxt != "=" and prev not in "=<>!":
                if stack and stack[-1]["at"] == 0 and ident:
                    stack[-1]["args"].append(ident)
            ident = ""
            i += 1
            continue
        # 空白字符不能清空 ident：`quiet = TRUE` 里 ident 与 `=` 之间隔着空格，
        # 早先一版在这里把 ident 清掉了，导致反例 gc(quiet = TRUE) 漏检。
        if not c.isspace():
            ident = ""
        i += 1
    return calls


def check_call_args(src):
    """返回所有「参数名不在该函数合法形参里」的调用说明。"""
    bad = []
    for line, name, args in _scan_calls(src):
        allowed = SIGNATURES[name]
        for a in args:
            if a not in allowed:
                bad.append(
                    "第 %d 行：%s() 没有名为 '%s' 的参数（合法形参：%s）"
                    % (line, name, a, ", ".join(sorted(allowed))))
    return bad


def check_file(path):
    """返回 (ok, info_dict)。ok 为 False 表示发现可疑问题。"""
    with open(path, "rb") as fh:
        raw = fh.read()
    src = raw.decode("utf-8")

    depth = {"(": 0, "{": 0, "[": 0}
    closer = {")": "(", "}": "{", "]": "["}
    errors = []
    i, n, line = 0, len(src), 1
    in_str = esc = in_comment = False

    while i < n:
        c = src[i]
        if c == "\n":
            line += 1
            in_comment = False
            i += 1
            continue
        if in_comment:
            i += 1
            continue
        if in_str:
            if esc:
                esc = False
            elif c == BACKSLASH:
                esc = True
            elif c == '"':
                in_str = False
            i += 1
            continue
        if c == "#":
            in_comment = True
            i += 1
            continue
        if c == '"':
            in_str = True
            i += 1
            continue
        if c in depth:
            depth[c] += 1
        elif c in closer:
            depth[closer[c]] -= 1
            if depth[closer[c]] < 0:
                errors.append("第 %d 行：多出一个 '%s'" % (line, c))
                depth[closer[c]] = 0
        i += 1

    # 扫一遍常见陷阱
    warnings = []
    for num, ln in enumerate(src.split("\n"), 1):
        s = ln.strip()
        # 只在「像 shell 命令」时才提示：cd 后面直接跟两个以上参数
        if re.match(r"^cd\s+\S+\s+\S+", s) and "<-" not in s:
            warnings.append("第 %d 行：cd 带了多个参数" % num)
        if "stopifnot" in s and "resultsNames" in s:
            warnings.append(
                "第 %d 行：用系数名做断言（因子水平一变就会误拦，应改为按数据做方向自检）"
                % num)

    crlf = raw.count(b"\r\n")
    has_bom = raw[:3] == b"\xef\xbb\xbf"
    calls = _scan_calls(src)
    bad_args = check_call_args(src)
    info = {
        "bytes": len(raw),
        "lines": src.count("\n") + 1,
        "paren": depth["("],
        "brace": depth["{"],
        "bracket": depth["["],
        "unterminated_string": in_str,
        "crlf": crlf,
        "bom": has_bom,
        "errors": errors,
        "warnings": warnings,
        "bad_args": bad_args,
        "calls_checked": len(calls),
        "sig_count": len({name for _, name, _ in calls}),
    }
    ok = (depth["("] == 0 and depth["{"] == 0 and depth["["] == 0
          and not in_str and not errors and crlf == 0 and not bad_args)
    return ok, info


def main(argv):
    if len(argv) < 2:
        print("用法: python tests/check_r_syntax.py <脚本路径>")
        return 2
    path = argv[1]
    if not os.path.exists(path):
        print("文件不存在:", path)
        return 2
    ok, info = check_file(path)
    print("文件        :", path)
    print("字节 / 行数 : %d / %d" % (info["bytes"], info["lines"]))
    print("圆括号 ( )  :", info["paren"])
    print("花括号 { }  :", info["brace"])
    print("方括号 [ ]  :", info["bracket"])
    print("字符串未闭合:", info["unterminated_string"])
    print("CRLF 行数   :", info["crlf"])
    print("含 BOM      :", info["bom"])
    print("已核对调用  : %d 处，涉及 %d 个已登记函数"
          % (info["calls_checked"], info["sig_count"]))
    for e in info["errors"][:10]:
        print("  !! 错误:", e)
    for b in info["bad_args"][:10]:
        print("  !! 参数名错误:", b)
    for w in info["warnings"][:10]:
        print("  [提示]", w)
    print()
    print("结论:", "括号与字符串均平衡、换行为纯 LF、参数字段名合法，未见明显语法错误"
          if ok else "疑似有问题，请检查")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
