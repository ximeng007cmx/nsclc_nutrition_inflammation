# -*- coding: utf-8 -*-
"""tools/md2html.py — 把本项目 manuscript/ 下的 Markdown 渲染成同款单文件 HTML。

用法：
  python tools/md2html.py <a.md> [<b.md> ...]       # 渲染指定文件（输出同名 .html）
  python tools/md2html.py --all                     # 渲染所有「已有同名 .html 的 .md」（不会凭空造新页面）
  python tools/md2html.py --all --roots manuscript,. # 指定搜索根（默认就是这两个）

样式（CSS）逐字取自本项目既有 HTML（manuscript/00_outline.html），
因此新渲染的页面与旧页面视觉完全一致。

支持的语法：
  # ## ### ####     标题
  > ...             引用块（可内嵌表格/列表，按内层 Markdown 递归渲染）
  ---               分隔线
  | a | b |         表格（含 |---| 分隔行）
  - / 1.            无序、有序列表
  内联：**加粗** *斜体* `code` ~~删除线~~ [文字](url)
"""

import html
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CSS = """:root{
  --bg:#f5f6f8; --card:#fff; --fg:#1a1d21; --fg2:#5b6472; --line:#e4e7ec;
  --accent:#0e6bb5; --soft:#eef5fc; --ok:#0f7b4f; --warn:#a35a00; --chip:#e9f1fa;
}
@media (prefers-color-scheme: dark){
  :root{
    --bg:#131619; --card:#1b1f23; --fg:#e8ecf0; --fg2:#a3adb8; --line:#2b3137;
    --accent:#68b1ea; --soft:#18232e; --ok:#5dd39e; --warn:#e0a458; --chip:#1f2c39;
  }
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font-family:"Microsoft YaHei","PingFang SC","Segoe UI",-apple-system,sans-serif;
  line-height:1.85;font-size:15.5px;-webkit-font-smoothing:antialiased}
.wrap{max-width:1000px;margin:0 auto;padding:36px 22px 90px}
h1{font-size:24px;line-height:1.5;margin:0 0 6px}
h2{font-size:19px;margin:44px 0 16px;padding:6px 0 6px 13px;border-left:4px solid var(--accent)}
h3{font-size:16px;margin:30px 0 12px;color:var(--accent)}
h4{font-size:15px;margin:20px 0 8px}
p{margin:0 0 14px;text-align:justify}
hr{border:none;border-top:1px solid var(--line);margin:34px 0}
strong{font-weight:700;background:linear-gradient(transparent 62%,var(--soft) 62%)}
code{background:var(--chip);padding:1px 6px;border-radius:5px;font-size:13.5px;
  font-family:Consolas,"Courier New",monospace}
blockquote{margin:16px 0;padding:12px 18px;border-left:3px solid var(--accent);
  background:var(--soft);border-radius:0 8px 8px 0;color:var(--fg2);font-size:14.5px}
ul,ol{padding-left:24px;margin:0 0 14px}
li{margin:5px 0}
.tw{overflow-x:auto;margin:16px 0 22px;border:1px solid var(--line);border-radius:10px}
table{border-collapse:collapse;width:100%;font-size:14px;background:var(--card)}
th{background:var(--chip);text-align:left;padding:10px 12px;font-weight:700;
  white-space:nowrap;border-bottom:1px solid var(--line)}
td{padding:9px 12px;border-top:1px solid var(--line);vertical-align:top}
tbody tr:nth-child(even){background:var(--soft)}
pre{background:var(--chip);border:1px solid var(--line);border-radius:10px;
  padding:14px 16px;overflow-x:auto;margin:16px 0 22px}
pre code{background:none;padding:0;font-size:13.5px;line-height:1.7}
.hero{background:var(--card);border:1px solid var(--line);border-radius:14px;
  padding:24px 26px;margin-bottom:8px}
.hero p{color:var(--fg2);font-size:14.5px}"""

HEAD = ('<!DOCTYPE html>\n<html lang="zh-CN">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<title>%s</title>\n<style>\n%s\n</style>\n</head>\n<body>\n<div class="wrap">\n')
TAIL = "</div>\n</body>\n</html>\n"


def inline(text):
    """行内标记。先转义，再把 code 段用占位符保护起来，避免内部被后续规则改写。"""
    text = html.escape(text, quote=True)
    codes = []

    def keep(m):
        codes.append(m.group(1))
        return "\x00%d\x00" % (len(codes) - 1)

    text = re.sub(r"`([^`]+)`", keep, text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*([^*\n]+?)\*(?!\*)", r"<em>\1</em>", text)
    text = re.sub(r"~~(.+?)~~", r"<s>\1</s>", text)
    text = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", r'<a href="\2">\1</a>', text)
    for i, c in enumerate(codes):
        text = text.replace("\x00%d\x00" % i, "<code>%s</code>" % c)
    return text


def split_row(line):
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [c.strip() for c in line.split("|")]


def is_sep(line):
    return bool(re.match(r"^\s*\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)*\|?\s*$", line))


def render(lines):
    out, i, n = [], 0, len(lines)
    while i < n:
        ln = lines[i]
        s = ln.strip()
        if not s:
            i += 1
            continue
        # 围栏代码块
        m = re.match(r"^\s*(`{3,}|~{3,})", ln)
        if m:
            ch = m.group(1)[0]
            code = []
            i += 1
            close = re.compile(r"^\s*" + re.escape(ch) + r"{3,}\s*$")
            while i < n and not close.match(lines[i]):
                code.append(lines[i])
                i += 1
            i += 1
            out.append("<pre><code>%s</code></pre>"
                       % html.escape("\n".join(code).strip("\n"), quote=True))
            continue
        # 分隔线
        if re.match(r"^(-{3,}|\*{3,}|_{3,})$", s):
            out.append("<hr>")
            i += 1
            continue
        # 标题
        m = re.match(r"^(#{1,6})\s+(.*)$", s)
        if m:
            lvl = min(len(m.group(1)), 4)
            out.append("<h%d>%s</h%d>" % (lvl, inline(m.group(2).strip()), lvl))
            i += 1
            continue
        # 引用块：收集连续的以 '>' 开头的行（含空 '>'），剥掉一层后递归
        if s.startswith(">"):
            inner = []
            while i < n and (lines[i].strip().startswith(">") or not lines[i].strip()):
                if not lines[i].strip():
                    # 引用块内部的空行只在后面仍跟着 '>' 时才算内部
                    j = i + 1
                    while j < n and not lines[j].strip():
                        j += 1
                    if j < n and lines[j].strip().startswith(">"):
                        inner.append("")
                        i += 1
                        continue
                    break
                inner.append(re.sub(r"^\s*>\s?", "", lines[i]))
                i += 1
            out.append("<blockquote>\n%s\n</blockquote>" % render(inner))
            i += 1 if (i < n and not lines[i].strip()) else 0
            continue
        # 表格
        if "|" in s and i + 1 < n and is_sep(lines[i + 1]):
            head = split_row(s)
            i += 2
            rows = []
            while i < n and lines[i].strip().startswith("|"):
                rows.append(split_row(lines[i]))
                i += 1
            buf = ["<div class='tw'><table>", "<thead><tr>"
                   + "".join("<th>%s</th>" % inline(c) for c in head) + "</tr></thead><tbody>"]
            for r in rows:
                r = r + [""] * (len(head) - len(r))
                buf.append("<tr>" + "".join("<td>%s</td>" % inline(c) for c in r) + "</tr>")
            buf.append("</tbody></table></div>")
            out.append("\n".join(buf))
            continue
        # 有序列表
        if re.match(r"^\d+\.\s+", s):
            items = []
            while i < n and re.match(r"^\s*\d+\.\s+", lines[i]):
                items.append(re.sub(r"^\s*\d+\.\s+", "", lines[i]).strip())
                i += 1
            out.append("<ol>" + "".join("<li>%s</li>" % inline(t) for t in items) + "</ol>")
            continue
        # 无序列表
        if re.match(r"^[-*+]\s+", s):
            items = []
            while i < n and re.match(r"^\s*[-*+]\s+", lines[i]):
                items.append(re.sub(r"^\s*[-*+]\s+", "", lines[i]).strip())
                i += 1
            out.append("<ul>" + "".join("<li>%s</li>" % inline(t) for t in items) + "</ul>")
            continue
        # 段落：连续非空、非块级起始的行，合并为一段
        para = []
        while i < n:
            t = lines[i].strip()
            if not t:
                break
            if (re.match(r"^(#{1,6})\s+", t) or t.startswith(">") or t.startswith("|")
                    or re.match(r"^(-{3,}|\*{3,}|_{3,})$", t)
                    or re.match(r"^(`{3,}|~{3,})", t)
                    or re.match(r"^\d+\.\s+", t) or re.match(r"^[-*+]\s+", t)):
                break
            para.append(t)
            i += 1
        if para:
            out.append("<p>%s</p>" % inline(" ".join(para)))
        else:
            i += 1
    return "\n".join(out)


def convert(md_path, html_path):
    src = open(md_path, encoding="utf-8").read()
    lines = src.split("\n")
    title = os.path.splitext(os.path.basename(md_path))[0]
    for ln in lines:
        m = re.match(r"^#\s+(.*)$", ln.strip())
        if m:
            title = m.group(1).strip()
            break
    body = render(lines)
    open(html_path, "w", encoding="utf-8").write(HEAD % (html.escape(title), CSS) + body + "\n" + TAIL)
    return title, len(body)


def main(argv):
    if not argv:
        print(__doc__)
        return 0
    roots = [os.path.join(ROOT, "manuscript"), ROOT]
    explicit = []
    if argv[0] == "--all":
        k = 1
        if k < len(argv) and argv[k] == "--roots":
            roots = [os.path.join(ROOT, p) for p in argv[k + 1].split(",")]
        for r in roots:
            if not os.path.isdir(r):
                continue
            for f in sorted(os.listdir(r)):
                if f.endswith(".md") and os.path.isfile(os.path.join(r, f[:-3] + ".html")):
                    explicit.append(os.path.join(r, f))
    else:
        explicit = [a if os.path.isabs(a) else os.path.join(ROOT, a) for a in argv]

    rc = 0
    for md in explicit:
        out = os.path.splitext(md)[0] + ".html"
        if not os.path.isfile(md):
            print("[SKIP] 不存在：%s" % md)
            rc = 1
            continue
        title, nb = convert(md, out)
        print("  %-46s → %s  （%d B 正文）" % (os.path.relpath(md, ROOT).replace("\\", "/"),
                                              os.path.basename(out), nb))
    print("共渲染 %d 个文件。" % len(explicit))
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
