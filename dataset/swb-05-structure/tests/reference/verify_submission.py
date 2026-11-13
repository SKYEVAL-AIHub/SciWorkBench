#!/usr/bin/env python
"""核验一份"根据 NMR/NOESY/MS 谱图数据推断有机分子结构"的提交件并判分。

用法：
    python answer/verify_submission.py <task_folder>

唯一需要读取的提交件是 `output/structures.csv`（三列：compound_id, smiles, inchikey，
compound_id 取 compound_A ~ compound_C）。两份评分表全部由本脚本机械判定，
不需要人工阅读自由文本，也不需要 LLM 参与判分：

- grading_no_tool.md（100 分）：分子式（15）+ 环数目（15）+ 官能团种类与数目（15）
  + 构造结构正确·不考虑立体（30）+ 立体结构正确·仅 compound_C（25）；
  一律由提交的 smiles 列重算后与参考值比较。
- grading_tool.md（100 分）：判据一致但把立体项压缩到 19 分，额外检查 inchikey 是否
  真由考生自己的 smiles 算出（2 分/化合物，共 6 分）。标准 InChIKey 无法手工推算，
  即"是否真的调用了工具"的证据。

立体考察范围：**只有 compound_C 考察立体构型**。判分口径：提交结构与参考结构为同一
立体异构体或互为对映异构体（全部手性中心取反；双键 E/Z 在镜像下保持不变）才算正确；
任何非对映异构体、部分指定或未指定立体的结构判错。compound_A/B 只考构造，
其提交 SMILES 写或不写立体均不影响判分。

退出码：0 = 无 FAIL，1 = 存在 FAIL，2 = 环境或参考文件缺失。

依赖：rdkit。
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors
except ImportError as exc:  # pragma: no cover
    print(f"环境缺少依赖：{exc}", file=sys.stderr)
    sys.exit(2)

RDLogger.DisableLog("rdApp.*")

STRUCT_COLS = ["compound_id", "smiles", "inchikey"]
ALL_CIDS = ["compound_A", "compound_B", "compound_C"]

# 与 answer/build_reference.py 保持完全一致的官能团定义
FG_SMARTS = {
    "ester_COO": "[CX3](=[OX1])[OX2]",          # 羧酸酯 C(=O)-O
    "alkene_C=C": "[#6X3;!R]=[#6X3;!R]",        # 非芳香、非环的烯烃 C=C（每条双键计 1）
    "alkyne_C#C": "[#6X2]#[#6X2]",              # 炔键（端炔亦命中）
    "epoxide": "C1OC1",                          # 环氧（三元环醚）
    "hydroxyl_OH": "[OX2H]",                     # 羟基（醇）
}
FG_TYPES = list(FG_SMARTS)

# 分值（每化合物）：得分点 -> {compound_id: 分值}。未列出的化合物表示该得分点不适用。
# 立体项仅 compound_C（A/B 不考立体）；InChIKey 为 grading_tool 独有。
W_NO_TOOL = {
    "一、分子式正确": {c: 5.0 for c in ALL_CIDS},
    "二、环数目正确": {c: 5.0 for c in ALL_CIDS},
    "三、官能团种类和数目正确": {c: 5.0 for c in ALL_CIDS},
    "四、构造结构正确（不考虑立体异构）": {c: 10.0 for c in ALL_CIDS},
    "五、立体结构正确": {"compound_C": 25.0},
}
W_TOOL = {
    "一、分子式正确": {c: 5.0 for c in ALL_CIDS},
    "二、环数目正确": {c: 5.0 for c in ALL_CIDS},
    "三、官能团种类和数目正确": {c: 5.0 for c in ALL_CIDS},
    "四、构造结构正确（不考虑立体异构）": {c: 10.0 for c in ALL_CIDS},
    "五、立体结构正确": {"compound_C": 19.0},
    "六、InChIKey 与自身结构一致": {c: 2.0 for c in ALL_CIDS},
}
# 得分点名称 -> compare() 结果键
ITEM_KEYS = {
    "一、分子式正确": "formula",
    "二、环数目正确": "rings",
    "三、官能团种类和数目正确": "fg",
    "四、构造结构正确（不考虑立体异构）": "constitution",
    "五、立体结构正确": "stereo",
    "六、InChIKey 与自身结构一致": "inchikey",
}

log: list[tuple[str, str]] = []


def rec(level: str, msg: str) -> None:
    log.append((level, msg))
    print(f"[{level}] {msg}")


# --------------------------------------------------------------------------- 化学

def parse(smiles: str):
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def fg_counts(mol) -> dict:
    return {name: len(mol.GetSubstructMatches(Chem.MolFromSmarts(s), uniquify=True))
            for name, s in FG_SMARTS.items()}


def stereo_element_count(mol) -> int:
    """已指定的立体元素数：已指定 E/Z 的双键 + 四面体手性中心（不含未指定）。

    注意顺序：必须先数 E/Z 再调用 FindMolChiralCenters
    （非 legacy 实现会清除双键立体标注）。
    """
    n = sum(1 for b in mol.GetBonds()
            if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ))
    n += sum(1 for _, c in Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                                     useLegacyImplementation=False)
             if c != "?")
    return n


def recompute(mol) -> dict:
    flat = Chem.Mol(mol)
    Chem.RemoveStereochemistry(flat)
    return {
        "formula": rdMolDescriptors.CalcMolFormula(flat),
        # 环数目按图论环数（键数-原子数+1），对桥环体系无 SSSR 歧义
        "rings": flat.GetNumBonds() - flat.GetNumAtoms() + 1,
        "fg": fg_counts(flat),
        # 构造基准：去立体后的标准 InChI（SMILES 字符串不是规范形式，InChI 才是）
        "constitution_inchi": Chem.MolToInchi(flat),
        # 立体基准：按提交 SMILES 原样（含其写的立体标注）的含立体标准 InChI
        "stereo_inchi": Chem.MolToInchi(mol),
        # 自洽性检查按提交 SMILES 原样重算（含考生自己写的立体标注），
        # 保证"用自己结构算出的 InChIKey"必然对得上
        "inchikey": Chem.MolToInchiKey(mol),
        "n_stereo_elements": stereo_element_count(mol),
    }


def compare(ref_entry: dict, mol) -> dict:
    """各得分点的布尔判定，全部由提交的 SMILES 重算后与参考值比较。

    立体项仅对参考条目标记 grades_stereo 的化合物（compound_C）计算；
    其余化合物 stereo 为 None（不适用，不参与判分）。
    """
    if mol is None:
        base = {"formula": False, "rings": False, "fg": False, "constitution": False}
        base["stereo"] = False if ref_entry["grades_stereo"] else None
        return base
    again = recompute(mol)
    ref_fg = ref_entry["functional_group_counts"]
    fg_detail = {t: (again["fg"][t], ref_fg[t]) for t in FG_TYPES}
    if ref_entry["grades_stereo"]:
        # 立体正确：必须指定参考结构的全部立体元素（防"不写立体碰运气"），
        # 且含立体 InChI 命中参考自身或其对映体之一
        stereo_ok = (
            again["n_stereo_elements"] >= ref_entry["n_stereo_elements"]
            and again["stereo_inchi"] in (ref_entry["stereo_standard_inchi"],
                                          ref_entry["enantiomer_standard_inchi"])
        )
    else:
        stereo_ok = None   # A/B 不考立体
    return {
        "formula": again["formula"] == ref_entry["molecular_formula"],
        "rings": again["rings"] == ref_entry["ring_count"],
        # 官能团"种类和数目"：被判的全部五类计数逐一相等才算对
        "fg": all(again["fg"][t] == ref_fg[t] for t in FG_TYPES),
        "constitution": again["constitution_inchi"] == ref_entry["constitution_standard_inchi"],
        "stereo": stereo_ok,
        "_recomputed": again,
        "_fg_detail": fg_detail,
    }


def norm_cid(raw: str) -> str:
    """宽容化 compound_id：去空白、统一大小写与分隔符。"""
    s = raw.strip().lower().replace(" ", "_").replace("-", "_")
    if s in {"a", "b", "c"}:
        return "compound_" + s.upper()
    if len(s) == 10 and s.startswith("compound_") and s[-1] in "abc":
        return "compound_" + s[-1].upper()
    return raw.strip()


# --------------------------------------------------------------------------- 主流程

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="核验 NMR/NOESY/MS 结构解析题提交件（全脚本判分）")
    ap.add_argument("task_folder", help="题目文件夹路径（其下应有 output/）")
    ap.add_argument("--reward-file", type=argparse.FileType("w"))
    args = ap.parse_args(argv)

    root = os.path.abspath(args.task_folder)
    struct_path = os.path.join(root, "output", "structures.csv")
    here = os.path.dirname(os.path.abspath(__file__))
    ref_path = os.path.join(here, "reference_values.json")
    if not os.path.isfile(ref_path):
        print(f"缺少参考文件：{ref_path}", file=sys.stderr)
        return 2
    reference = json.load(open(ref_path, encoding="utf-8"))

    print("=" * 78)
    print("第 1 部分  提交件读取（不单独计分；缺失或不合规会使相应得分点为 0）")
    print("=" * 78)

    submitted = {}
    if not os.path.isfile(struct_path):
        rec("FAIL", "缺少 output/structures.csv，两份评分表的全部得分点均为 0")
    else:
        with open(struct_path, newline="", encoding="utf-8-sig") as f:
            rows = [r for r in csv.reader(f) if any(c.strip() for c in r)]
        header = [c.strip() for c in rows[0]] if rows else []
        if header == STRUCT_COLS:
            rec("PASS", "structures.csv 表头与规定的 3 列一致")
        else:
            rec("WARN", f"structures.csv 表头为 {header}，与规定的 {STRUCT_COLS} 不一致；"
                        "仍按其中可识别的列判分")
        for r in rows[1:]:
            d = dict(zip(header, [c.strip() for c in r]))
            cid = norm_cid(d.get("compound_id", ""))
            if cid in ALL_CIDS:
                if cid in submitted:
                    rec("WARN", f"{cid} 出现多行，按第一行判分")
                else:
                    submitted[cid] = d
            elif d.get("compound_id"):
                rec("WARN", f"无法识别的 compound_id：{d['compound_id']!r}，该行不计分")
        missing = [c for c in ALL_CIDS if c not in submitted]
        if missing:
            rec("FAIL", f"structures.csv 缺少以下 compound_id 的行：{missing}")
        else:
            rec("PASS", "structures.csv 含三个 compound_id 共 3 行")

    mols, bad = {}, []
    for cid in ALL_CIDS:
        mols[cid] = parse(submitted.get(cid, {}).get("smiles", ""))
        if mols[cid] is None:
            bad.append(cid)
    if bad:
        rec("FAIL", f"以下化合物的 SMILES 缺失或无法解析，其得分点为 0：{bad}")
    else:
        rec("PASS", "3 个 SMILES 全部可正常解析")

    # ---------------------------------------------------------- 得分点一~五：三个化合物
    print()
    print("=" * 78)
    print("第 2 部分  得分点一~五 · 三个化合物的结构比对（按 compound_id 直接对应；"
          "立体项仅 compound_C）")
    print("=" * 78)

    flags = {}
    for cid in ALL_CIDS:
        entry = reference["compounds"][cid]
        mol = mols.get(cid)
        cmp = compare(entry, mol)
        flags[cid] = cmp
        if mol is None:
            rec("WARN", f"{cid}: 结构缺失或无法解析，全部得分点不给分")
            continue
        again = cmp["_recomputed"]
        detail = (f"分子式 {again['formula']}/{entry['molecular_formula']}、"
                  f"环数 {again['rings']}/{entry['ring_count']}")
        marks = []
        for name, key in list(ITEM_KEYS.items())[:5]:
            v = cmp[key]
            if v is None:
                marks.append(f"–{name.split('、')[1]}(不适用)")
            else:
                marks.append(f"{'✓' if v else '✗'}{name.split('、')[1]}")
        if entry["grades_stereo"]:
            detail += (f"、立体元素 {again['n_stereo_elements']}"
                       f"/{entry['n_stereo_elements']}")
        rec("PASS" if all(v for v in cmp.values()
                          if isinstance(v, bool) and v is not None)
            and all(cmp[k] for k in ("formula", "rings", "fg", "constitution"))
            else "WARN",
            f"{cid}（文献 {entry['source_label']}）: " + " ".join(marks) + f"（{detail}）")
        if not cmp["fg"]:
            fg_note = "、".join(f"{t} {got}/{want}"
                                for t, (got, want) in cmp["_fg_detail"].items() if got != want)
            print(f"        · 官能团计数不一致项：{fg_note}")
        if not cmp["constitution"]:
            print(f"        · 去立体 InChI：提交 {again['constitution_inchi']}")
            print(f"        · 去立体 InChI：参考 {entry['constitution_standard_inchi']}")
        if entry["grades_stereo"] and cmp["constitution"] and not cmp["stereo"]:
            if again["n_stereo_elements"] < entry["n_stereo_elements"]:
                print(f"        · 立体判错原因：未指定全部 {entry['n_stereo_elements']} 个立体元素")
            else:
                print("        · 立体判错原因：与参考既不同一也非对映体关系（非对映异构体或部分立体错误）")
            print(f"        · 含立体 InChI：提交 {again['stereo_inchi']}")
            print(f"        · 含立体 InChI：参考 {entry['stereo_standard_inchi']}")

    # ---------------------------------------------------------- 得分点六：InChIKey 自洽性
    print()
    print("=" * 78)
    print("第 3 部分  InChIKey 自洽性（grading_tool.md 独有的得分点）")
    print("=" * 78)
    ik_ok = {}
    for cid in ALL_CIDS:
        mol = mols.get(cid)
        if mol is None:
            ik_ok[cid] = False
            rec("WARN", f"{cid}: 结构缺失或无法解析，无法重算 InChIKey，本项不给分")
            continue
        want = recompute(mol)["inchikey"]
        raw = submitted.get(cid, {}).get("inchikey", "").strip()
        ik_ok[cid] = raw.upper() == want.upper()
        rec("PASS" if ik_ok[cid] else "WARN",
            f"{cid}: 申报 InChIKey {raw or '（空）'}，由提交 SMILES 重算得 {want}"
            + ("" if ik_ok[cid] else "  → 不一致"))
    if sum(ik_ok.values()) == 0:
        rec("FAIL", "没有任何一个化合物的申报 InChIKey 与其自身 SMILES 相符，"
                    "无法认定考生调用过化学信息学工具")
    for cid in ALL_CIDS:
        flags[cid]["inchikey"] = ik_ok[cid]

    # --------------------------------------------------------------- 判分汇总
    def build(weights: dict):
        table = []
        for name, wmap in weights.items():
            key = ITEM_KEYS[name]
            full = sum(wmap.values())
            got = sum(w for c, w in wmap.items() if flags[c].get(key))
            hit = [c for c in wmap if flags[c].get(key)]
            table.append((name, full, got, f"{len(hit)}/{len(wmap)} 化合物：{hit}"))
        return table

    for title, weights in (("grading_no_tool.md", W_NO_TOOL), ("grading_tool.md", W_TOOL)):
        print()
        print("=" * 78)
        print(f"判分结果（{title}，满分 100）")
        print("=" * 78)
        total = 0.0
        for name, full, got, note in build(weights):
            total += got
            print(f"  {name:<36} {got:>6.2f} / {full:<4g}  {note}")
        print(f"  {'合计':<36} {total:>6.2f} / 100")
        if title == "grading_no_tool.md" and args.reward_file:
            args.reward_file.write(f"{total / 100:.8f}\n")

    print()
    print("提示：以上两份判分结果全部由脚本得出。得分点一~五一律以提交的 smiles 列重算为准，"
          "不采用任何申报值；立体项仅对 compound_C 判分，且当且仅当与参考为同一立体异构体"
          "或互为对映异构体；compound_A/B 写或不写立体均不影响判分；"
          "InChIKey 得分点考察申报值是否真由该 smiles 算出。")
    fails = [m for lv, m in log if lv == "FAIL"]
    print(f"FAIL {len(fails)} 项，WARN {len([1 for lv, _ in log if lv == 'WARN'])} 项")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
