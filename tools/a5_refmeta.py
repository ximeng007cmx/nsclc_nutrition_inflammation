# -*- coding: utf-8 -*-
"""A5 引用台账：按拟引 PMID 拉取 Europe PMC 全元数据 + CrossRef 双源核验。
只读远端，写出 manuscript/refs/a5_refmeta.csv / .json。
铁律：任何一条若 EPMC 查不到 / DOI 在 CrossRef 不存在，就地报出，绝不静默。
"""
import csv, json, os, sys, time, urllib.request, urllib.parse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, "manuscript", "refs")
UA = {"User-Agent": "ali-pathA-refmeta/1.0 (mailto:chenmx@example.org)"}

# 角色分组（决定引用位置；不改动任何数值）
SLATE = [
    # A. ALI/PNI/NLR 在 NSCLC 的临床预后证据（临床侧锚点）
    ("23530866", "A", "ALI 在 NSCLC 的原始预后报道（Jafri 2013 BMC Cancer）"),
    ("32010567", "A", "NLR/LMR/PLR/ALI 在 IV 期 NSCLC 的预后比较"),
    ("35833264", "A", "ALI 为肺癌 OS 最优炎症标志物"),
    ("31677642", "A", "ALI 在肿瘤患者预后中的 meta 分析"),
    ("35297277", "A", "ALI 在 NSCLC 的 meta 分析"),
    ("29049691", "A", "ALI 是 NSCLC 术后独立预后因素"),
    ("36835030", "A", "炎症/营养评分「预后但非预测」——支持我们「不是预测工具」的定位"),
    ("37186387", "A", "全身炎症与根治性切除 NSCLC 预后相关并伴随免疫抑制微环境（最贴近的临床-免疫桥）"),
    ("34481329", "A", "ALI 与 NSCLC 免疫检查点抑制剂疗效关联"),
    ("36980740", "A", "可切除 NSCLC 中炎症指数的预后价值"),
    # B. 营养-炎症的分子/生理学基础（为什么白蛋白不是「营养」而是负性急性期蛋白）
    ("15660573", "B", "白蛋白：炎症与营养的关系（经典综述）"),
    ("34181998", "B", "血清白蛋白与炎症"),
    ("3584463",  "B", "炎症时肝合成白蛋白下降的分子机制（经典）"),
    ("35909571", "B", "脂肪因子/肝因子/肌因子与脂肪组织炎症"),
    ("33846641", "B", "秋水仙碱经肝诱导肝因子抑制髓系细胞活化（Nat Metab）"),
    ("33436955", "B", "LECT2 作为肝因子经组织巨噬细胞连接脂肪肝与炎症（Sci Rep）"),
    ("26699087", "B", "富组氨酸糖蛋白促进慢性肝病巨噬细胞活化与炎症"),
    ("38344945", "B", "NLR 与营养参数、瘤内免疫谱及 pembrolizumab 结局相关（头颈鳞癌）"),
    # C. 最接近的前人工作
    ("25238252", "C", "全身炎症 + 营养状态 + 肿瘤免疫微环境共同决定切除 NSCLC 结局（最接近的前人工作）"),
    # D. 表达签名可重复性 / 过拟合 / 方法学
    ("15705458", "D", "微阵列预测结局：多重随机验证策略（Michiels 2005 Lancet）"),
    ("15308542", "D", "结局签名基因是否存在唯一集合（Ein-Dor 2005 Bioinformatics）"),
    ("22028643", "D", "多数随机表达签名与乳腺癌结局显著相关（Venet 2011 PLoS Comput Biol）"),
    ("8668867",  "D", "多变量预后模型：开发、假设检验与误差（Harrell 1996 Stat Med）"),
    # E. 肿瘤纯度混杂
    ("24113773", "E", "由表达谱推断肿瘤纯度与间质/免疫细胞混合（Yoshihara 2013 ESTIMATE，Nat Commun）"),
    ("26634437", "E", "泛癌肿瘤纯度系统分析（Aran 2015 Nat Commun）——「肿瘤纯度混杂」的正主"),
    ("29141981", "E", "肿瘤纯度对跨癌种免疫基因表达与聚类的影响（Cancer Immunol Res）"),
    # F. 瘤内 T 细胞浸润的预后背景
    ("16421594", "F", "CD8+ 与 CD4+ T 细胞共浸润为 NSCLC 的有利预后因素（Br J Cancer，经典）"),
    ("26845192", "F", "肿瘤浸润淋巴细胞在 NSCLC 发生、进展与预后中的作用"),
    # G. 炎症表型与免疫治疗
    ("33220234", "G", "拮抗性炎症表型决定肿瘤命运与 ICB 应答（Immunity 2020）"),
    ("31684954", "G", "肿瘤炎症签名（TIS）与抗 PD-1 获益相关（CERTIM 泛癌队列）"),
    ("31683225", "G", "肿瘤炎症与 EMT 基因签名高精度预测肺癌 ICB 应答"),
    ("35398717", "G", "全身炎症标志物预测 NSCLC 一线 ICI 生存"),
    # H. TCGA 单队列「签名洪流」语境（T4/T6 的代言，仅少量）
    ("37657358", "H", "肺腺癌二硫死亡相关基因预后签名（TCGA 单队列类型的代表）"),
    ("34844602", "H", "NSCLC 免疫分子亚型预后刻画（TCGA 单队列类型的代表）"),
    # I. TCGA / Toil / 方法学复现硬需求
    ("25079552", "I", "TCGA-LUAD 分子分型（Nature 2014）"),
    ("22960745", "I", "TCGA-LUSC 基因组刻画（Nature 2012）"),
    ("29625055", "I", "TCGA 泛癌临床数据资源（Liu 2018 Cell）"),
    ("28398314", "I", "Toil 可重复的开放生物医学数据分析（Vivian 2017 Nat Biotechnol）"),
    ("26209429", "I", "TCGA RNA-Seq 替代预处理改进分析结果（Rahman 2015 Bioinformatics）"),
    ("32444850", "I", "Xena 平台（Goldman 2020 Nat Biotechnol）"),
    ("19847166", "I", "ssGSEA 原始文献（Barbie 2009 Nature）"),
    ("16199517", "I", "GSEA 原始文献（Subramanian 2005 PNAS）"),
    ("23323831", "I", "GSVA（Hänzelmann 2013 BMC Bioinformatics）"),
    ("26771021", "I", "MSigDB hallmark 基因集（Liberzon 2015 Cell Syst）"),
    ("20808728", "I", "glmnet 正则化路径（Friedman 2010 J Stat Softw）"),
    ("BENJAMINI1995", "I", "Benjamini-Hochberg FDR 方法（1995 JRSS-B，PubMed 未收录）"),
]

DOI_OVERRIDE = {
    "BENJAMINI1995": "10.1111/j.2517-6161.1995.tb02031.x",
    "20808728": "10.18637/jss.v033.i01",
}


def epmc_core(pmid):
    u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
        {"query": "EXT_ID:%s AND SRC:MED" % pmid, "format": "json",
         "resultType": "core", "pageSize": 1})
    with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=45) as r:
        d = json.load(r)
    if not d.get("hitCount"):
        return None
    return d["resultList"]["result"][0]


def crossref(doi):
    u = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=45) as r:
            m = json.load(r)["message"]
        return True, (m.get("container-title") or [""])[0], (m.get("title") or [""])[0]
    except Exception as e:
        return False, None, str(e)[:120]


def main():
    rows, problems = [], []
    total = len(SLATE)
    for idx, (pmid, grp, role) in enumerate(SLATE, 1):
        print("  [%2d/%d] %s" % (idx, total, pmid), flush=True)
        if pmid == "BENJAMINI1995":
            doi = DOI_OVERRIDE[pmid]
            ok, jr, tt = crossref(doi)
            rows.append(dict(pmid="(none)", doi=doi, group=grp, role=role,
                             authors="Benjamini Y, Hochberg Y.", title=tt or "",
                             journal="J R Stat Soc Series B Stat Methodol", year="1995",
                             vol="57", iss="1", pages="289-300",
                             epmc="not_indexed", crossref="ok" if ok else "FAIL"))
            if not ok:
                problems.append((pmid, "CrossRef 查不到 %s" % doi))
            time.sleep(0.3)
            continue

        r = epmc_core(pmid)
        if r is None:
            problems.append((pmid, "Europe PMC 查不到该 PMID"))
            continue
        ji = r.get("journalInfo", {}) or {}
        j = ji.get("journal", {}) or {}
        doi = (r.get("doi") or "").strip() or DOI_OVERRIDE.get(pmid, "")
        cr_ok, cr_jr, cr_tt = crossref(doi) if doi else (None, None, "no doi")
        if cr_ok is False:
            problems.append((pmid, "CrossRef 查不到 DOI %s" % doi))
        rows.append(dict(
            pmid=pmid, doi=doi, group=grp, role=role,
            authors=(r.get("authorString") or "").strip(),
            title=(r.get("title") or "").strip().rstrip("."),
            journal=j.get("title") or "", year=r.get("pubYear") or "",
            vol=ji.get("volume") or "", iss=ji.get("issue") or "",
            pages=(r.get("pageInfo") or ""),
            epmc="ok", crossref=("ok" if cr_ok else ("FAIL" if cr_ok is False else "n/a")),
        ))
        time.sleep(0.3)

    os.makedirs(OUTDIR, exist_ok=True)
    csvp = os.path.join(OUTDIR, "a5_refmeta.csv")
    with open(csvp, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    jsonp = os.path.join(OUTDIR, "a5_refmeta.json")
    with open(jsonp, "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=1)

    print("台账 %d 条（含 1 条 PubMed 未收录的 B-H）" % len(rows))
    for r in rows:
        print("%-11s %-4s %-30s %s" % (r["pmid"], r["year"],
                                       (r["journal"] or "")[:30], (r["title"] or "")[:70]))
    print()
    if problems:
        print("⚠ 问题 %d 条：" % len(problems))
        for p, m in problems:
            print("   %s  %s" % (p, m))
    else:
        print("✅ 全部通过：EPMC 归属 OK，CrossRef DOI 全部存在")
    print("\n→ %s\n→ %s" % (csvp, jsonp))


if __name__ == "__main__":
    main()
