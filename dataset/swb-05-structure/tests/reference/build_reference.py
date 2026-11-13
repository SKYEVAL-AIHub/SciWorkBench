#!/usr/bin/env python
"""生成本题的参考答案产物（可复现）。

用法：
    python answer/build_reference.py

输出（写入本脚本所在目录）：
    reference_values.json      机器可读参考值与判分参数（含参考 SMILES、去立体/含立体
                               标准 InChI、对映体 InChI、分子式、环数、官能团计数等
                               全部得分点的参考值）

题目：三个有机化合物（编号 compound_A ~ compound_C，对应文献中 12a / 12b / 13，
见 reference_report.md）。A/B 为带环氧、端炔与 (E)-烯酸甲酯侧链的多取代脂环醇；
C 是 A 经 Cp2TiCl 自由基串环化得到的三环产物（环外亚甲基 + 饱和甲酯，比 12b 的
串环化产物少一个桥亚甲基），
其相对构型需由随数据给出的 NOESY 谱（data/compound_C_noesy.png）确定。

**立体考察范围**：只有 compound_C 考察立体构型（参考 SMILES 及其对映异构体算正确，
非对映异构体算错误）；compound_A/B 只考构造，其谱图不含足够构型信息，立体不参与判分。

判分全部从提交的 smiles 列重算后与本文件生成的参考值比较（本题 3 个化合物，分值按满分 100）：
    一、分子式（5 分/化合物）；二、环数目（5 分/化合物）；
    三、官能团种类与数目（5 分/化合物）；四、构造结构正确（不考虑立体异构，10 分/化合物）；
    五、立体结构正确（仅 compound_C，no_tool 25 / tool 19 分）；
    六、InChIKey 自洽性（tool 独有，2 分/化合物）。

依赖：rdkit。
"""

from __future__ import annotations

import json
import os

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))

PROTON = 1.007276    # Da

# compound_id -> (文献编号, 数据文件, 附加数据文件, HRMS 实测 m/z, 离子类型, 是否考立体, 参考 SMILES)
REFERENCE = [
    ("compound_A", "12a", "data/compound_A.csv", [], 264.1362, "[M]+", False,
     r"O[C@@H]1[C@@](C)(C/C=C/C(OC)=O)C[C@@](O2)(CC#C)[C@@]2([H])C1"),
    ("compound_B", "12b", "data/compound_B.csv", [], 278.1511, "[M]+", False,
     r"O[C@@H]1[C@@](C)(CC/C=C/C(OC)=O)C[C@@](O2)(CC#C)[C@@]2([H])C1"),
    ("compound_C", "13", "data/compound_C.csv", ["data/compound_C_noesy.png"],
     267.1597, "[M+H]+", True,
     r"O[C@@H]1[C@@](C)(C[C@]2([H])[C@@H](C(OC)=O)C(C3)=C)C[C@]23[C@H](O)C1"),
]

# 被判的官能团类型：名称 -> SMARTS。计数一律用 uniquify 子结构匹配。
# alkene 只计非环 C=C（C 的环外亚甲基 =CH2 中 =CH2 碳非 X3，且稠环桥头双键
# 属环内双键，均不计入；A/B 的 (E)-侧链烯键计 1）。
FG_SMARTS = {
    "ester_COO": "[CX3](=[OX1])[OX2]",          # 羧酸酯 C(=O)-O
    "alkene_C=C": "[#6X3;!R]=[#6X3;!R]",        # 非芳香、非环的烯烃 C=C（每条双键计 1）
    "alkyne_C#C": "[#6X2]#[#6X2]",              # 炔键（端炔亦命中）
    "epoxide": "C1OC1",                          # 环氧（三元环醚）
    "hydroxyl_OH": "[OX2H]",                     # 羟基（醇）
}


def fg_counts(mol) -> dict:
    return {name: len(mol.GetSubstructMatches(Chem.MolFromSmarts(s), uniquify=True))
            for name, s in FG_SMARTS.items()}


def enantiomer_inchi(mol) -> str:
    """对 mol 的镜像（全部四面体手性中心取反）计算标准 InChI。

    双键 E/Z 在镜像操作下保持不变（E 的对映体仍是 E，Z 亦然），故不翻转。
    标准 InChI 的 /m 标志区分对映体：镜像的 /t 层相同、/m 取反。
    无任何手性中心时返回与自身相同的 InChI（此时"对映体"即自身）。
    """
    en = Chem.Mol(mol)
    for a in en.GetAtoms():
        if a.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                                Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            a.InvertChirality()
    Chem.AssignStereochemistry(en, cleanIt=True, force=True)
    return Chem.MolToInchi(en)


def main() -> None:
    data = {"compounds": {}, "compound_order": [], "fg_types": list(FG_SMARTS)}

    for cid, source, datafile, extra, mz, ion, stereo, smi in REFERENCE:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"参考 SMILES 无法解析：{smi}")
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        # 注意：必须先数 E/Z 再调用 FindMolChiralCenters（非 legacy 实现会清除双键立体标注）
        n_ez = sum(1 for b in mol.GetBonds()
                   if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ))
        stereo_inchi = Chem.MolToInchi(mol)
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                            useLegacyImplementation=False)
        n_centers = sum(1 for _, c in centers if c != "?")
        canonical = Chem.MolToSmiles(mol)
        flat = Chem.Mol(mol)
        Chem.RemoveStereochemistry(flat)
        # 构造判分基准：去立体后的标准 InChI（SMILES 字符串不是规范形式，InChI 才是）
        flat_inchi = Chem.MolToInchi(flat)
        data["compounds"][cid] = {
            "source_label": source,
            "data_file": datafile,
            "extra_data_files": extra,
            "hrms_found_mz": mz,
            "hrms_ion_type": ion,
            "grades_stereo": stereo,
            "reference_smiles": smi,
            "canonical_smiles": canonical,
            # 得分点四（构造，不考虑立体）基准
            "constitution_standard_inchi": flat_inchi,
            # 得分点五（立体，仅 compound_C）基准：参考自身与对映体的含立体标准 InChI
            "stereo_standard_inchi": stereo_inchi,
            "enantiomer_standard_inchi": enantiomer_inchi(mol),
            "molecular_formula": rdMolDescriptors.CalcMolFormula(flat),
            "monoisotopic_mass_Da": round(Descriptors.ExactMolWt(flat), 4),
            "average_molar_mass_g_per_mol": round(Descriptors.MolWt(flat), 2),
            "heavy_atom_count": flat.GetNumHeavyAtoms(),
            # 环数目按图论环数（键数-原子数+1），对桥环体系无 SSSR 歧义
            "ring_count": flat.GetNumBonds() - flat.GetNumAtoms() + 1,
            "functional_group_counts": fg_counts(flat),
            "n_stereocenters_specified": n_centers,
            "n_double_bond_stereo_specified": n_ez,
            "n_stereo_elements": n_centers + n_ez,
        }
        data["compound_order"].append(cid)

    data["grading_notes"] = {
        "output_file": "output/structures.csv",
        "output_columns": ["compound_id", "smiles", "inchikey"],
        "grading_no_tool": {
            "分子式正确": {"per_compound": 5, "source": "由提交 smiles 重算 CalcMolFormula 与参考比较"},
            "环数目正确": {"per_compound": 5,
                          "source": "图论环数（键数-原子数+1），参考值 A=2, B=2, C=3"},
            "官能团种类和数目正确": {"per_compound": 5,
                                    "source": "ester_COO/alkene_C=C/alkyne_C#C/epoxide/hydroxyl_OH 五类计数全部相等"},
            "构造结构正确（不考虑立体异构）": {"per_compound": 10,
                                             "source": "去立体后的标准 InChI 逐字比对"},
            "立体结构正确": {"compound_C": 25,
                            "source": "仅 compound_C；含立体标准 InChI 等于参考或其对映体（全部手性中心取反；"
                                      "E/Z 在镜像下不变）；非对映异构体、部分指定或未指定立体均判错"},
        },
        "grading_tool": {
            "分子式正确": {"per_compound": 5},
            "环数目正确": {"per_compound": 5},
            "官能团种类和数目正确": {"per_compound": 5},
            "构造结构正确（不考虑立体异构）": {"per_compound": 10},
            "立体结构正确": {"compound_C": 19},
            "InChIKey与自身结构一致": {"per_compound": 2},
        },
        "stereo_policy": "本题只有 compound_C 考察立体异构：立体结构正确当且仅当提交结构与参考结构"
                         "为同一立体异构体或互为对映异构体（外消旋合成无法区分对映体）；"
                         "任何非对映异构体判错。compound_C 的相对构型由随题给出的 NOESY 谱确定。"
                         "compound_A/B 的谱图不含足够构型信息，只考构造，提交中写或不写立体均不影响判分。",
        "hrms_policy": "三个化合物均给出 HRMS 实测 m/z 与离子类型（A/B 为 FI/FD 的 [M]+，"
                       "C 为 FI 的 [M+H]+），但刻意不给出文献中 calcd 的元素组成，"
                       "分子式需由精确质量结合 NMR 计数推出。",
    }

    with open(os.path.join(HERE, "reference_values.json"), "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    for cid in data["compound_order"]:
        c = data["compounds"][cid]
        adduct = PROTON if c["hrms_ion_type"] == "[M+H]+" else 0.0
        calc = c["monoisotopic_mass_Da"] + adduct
        dev = (calc - c["hrms_found_mz"]) / c["hrms_found_mz"] * 1e6
        print(f"{cid}: {c['molecular_formula']}  rings={c['ring_count']}  "
              f"grades_stereo={c['grades_stereo']}  "
              f"stereo={c['n_stereocenters_specified']} centers + {c['n_double_bond_stereo_specified']} E/Z")
        print(f"    官能团 {c['functional_group_counts']}")
        print(f"    HRMS {c['hrms_ion_type']} calc {calc:.4f} / found {c['hrms_found_mz']} → {dev:+.1f} ppm")


if __name__ == "__main__":
    main()
