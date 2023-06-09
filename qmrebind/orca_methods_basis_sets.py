from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfMerger
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# QM2 Methods
df1 = pd.DataFrame(["AM1", "PM3"])
df1.columns = ["QM2 Semiempirical Methods"]
df2 = pd.DataFrame(["XTB", "XTB1"])
df2.columns = ["Tight Binding DFT Methods"]
df3 = pd.DataFrame(["HF-3C", "PBEH-3C", "R2SCAN-3C"])
df3.columns = ["Composite Methods"]
df_qm2_methods = pd.concat([df1, df2, df3], axis=1)
df_qm2_methods.fillna("-", inplace=True)
df_qm2_methods.reset_index(drop=True, inplace=True)
# QM Basis Sets
df1 = pd.DataFrame(
    ["STO-3G", "3-21G", "3-21GSP", "4-22GSP", "6-31G", "m6-31G", "6-311G"]
)
df1.columns = ["Pople Basis Sets"]
df2 = pd.DataFrame(
    [
        "6-31G(d)",
        "6-31G(d,p)",
        "6-31G(2d)",
        "6-31G(2df)",
        "6-31G(2d,p)",
        "6-31G(2d,2p)",
        "6-31G(2df,2p)",
        "6-31G(2df,2pd)",
        "6-311G(d)",
        "6-311G(d,p)",
        "6-311G(2d)",
        "6-311G(2d,p)",
        "6-311G(2d,2p)",
        "6-311G(2df,2pd)",
        "6-311G(3df)",
        "6-311G(3df,3pd)",
    ]
)
df2.columns = ["Pople Polarized Basis Sets"]
df3 = pd.DataFrame(
    [
        "6-31+G(d)",
        "6-31+G(d,p,)",
        "6-31+G(2d)",
        "6-31+G(2df)",
        "6-31+G(2d,p)",
        "6-31+G(2d,2p)",
        "6-31+G(2df,2pd)",
        "6-311+G(d)",
        "6-311+G(d,p)",
        "6-311+G(2d)",
        "6-311+G(2df)",
        "6-311+G(2d,p)",
        "6-311+G(2d,2p)",
        "6-311+G(2df,2p)",
        "6-311+G(2df.2pd)",
        "6-311+G(3df)",
        "6-311+G(3df,3pd)",
        "6-31++G(d,p)",
        "6-31++G(2d,p)",
        "6-31++G(2d,2p)",
        "6-31++G(2df,2p)",
        "6-31++G(2df,2pd)",
        "6-311++G(d,p)",
        "6-311++G(2d,p)",
        "6-311++G(2d,2p)",
        "6-311++G(2df,2p)",
        "6-311++G(2df,2pd)",
        "6-311++G(3df,3pd)",
    ]
)
df3.columns = ["Pople Polarized Diffused Basis Sets"]
df4 = pd.DataFrame(
    [
        "cc-pVDZ",
        "cc-pVTZ",
        "cc-pVQZ",
        "cc-pV5Z",
        "cc-pV6Z",
        "aug-cc-pVDZ",
        "aug-cc-pVTZ",
        "aug-cc-pVDZ",
        "aug-cc-pV5Z",
        "aug-cc-pV6Z",
        "cc-pCVDZ",
        "cc-pCVTZ",
        "cc-pCVQZ",
        "cc-pCV5Z",
        "cc-pCV6Z",
        "aug-cc-pCVDZ",
        "aug-cc-pCVTZ",
        "aug-cc-pCVQZ",
        "aug-cc-pCV5Z",
        "aug-cc-pCV6Z",
        "cc-pwCVDZ",
        "cc-pwCVTZ",
        "cc-pwCVQZ",
        "cc-pwCV5Z",
        "aug-cc-pwCVDZ",
        "aug-cc-pwCVTZ",
        "aug-cc-pwCVQZ",
        "aug-cc-pwCV5Z",
        "cc-pVD(+d)Z",
        "cc-pVT(+d)Z",
        "cc-pVQ(+d)Z",
        "cc-pV5(+d)Z",
    ]
)
df4.columns = ["Correlation Consistent Basis Sets"]
df5 = pd.DataFrame(
    [
        "def2-SVP",
        "def2-SVP(P)",
        "def2-TZVP",
        "def2-TZVP(-f)",
        "def2-TZVPP",
        "def2-QZVP",
        "def2-QZVPP",
        "SV",
        "SV(P)",
        "SVP",
        "TZV",
        "TZV(P)",
        "TZVP",
        "TZVPP",
        "QZVP",
        "QZVPP",
    ]
)
df5.columns = ["DEF2 Basis Sets"]
df6 = pd.DataFrame(
    [
        "ma-def2-SVP",
        "ma-def2-SV(P)",
        "ma-def2-TZVP",
        "ma-def2-TZVP(-f)",
        "ma-def2-TZVPP",
        "ma-def2-QZVPP",
        "def2-SVPD",
        "def2-TZVPD",
        "def2-TZVPPD",
        "def2-QZVPD",
        "def2-QZVPPD",
    ]
)
df6.columns = ["DEF2 Diffused Basis Sets"]
df7 = pd.DataFrame(["Def2/J", "f2/JKsmall", "x2c/J"])
df7.columns = ["Auxiliary Coulomb Fit Basis Sets"]
df8 = pd.DataFrame(
    [
        "Def2/JK",
        "Def2/JKsmall",
        "cc-pVTZ/JK",
        "cc-pVQZ/JK",
        "cc-pV5Z/JK",
        "aug-cc-pVTZ/JK",
        "aug-cc-pV5Z/JK",
    ]
)
df8.columns = ["Auxiliary Coulomb Fit / Exchange Basis Sets"]
df9 = pd.DataFrame(
    [
        "Def2-SVP/C",
        "Def2-TZVP/C",
        "Def2-TZVPP/C",
        "Def2-QZVPP/C",
        "Def2-SVPD/C",
        "Def2-TZVPD/C",
        "Def2-TZVPPD/C",
        "Def2-QZVPPD/C",
        "cc-pVDZ/C",
        "cc-pVTZ/C",
        "cc-pVQZ/C",
        "cc-pV5Z/C",
        "cc-pV6Z/C",
        "aug-cc-pVDZ/C",
        "aug-cc-pVTZ/C",
        "aug-cc-pVQZ/C",
        "aug-cc-pV5Z/C",
        "aug-cc-pV6Z/C",
        "cc-pwCVDZ/C",
        "cc-pwCVTZ/C",
        "cc-pwCVQZ/C",
        "cc-pwCV5Z/C",
        "aug-cc-pwCVDZ/C",
        "aug-cc-pwCVTZ/C",
        "aug-cc-pwCVQZ/C",
        "aug-cc-pwCV5Z/C",
        "cc-pVDZ-PP/C",
        "cc-pVTZ-PP/C",
        "cc-pVQZ-PP/C",
        "aug-cc-pVDZ-PP/C",
        "aug-cc-pVTZ-PP/C",
        "aug-cc-pVQZ-PP/C",
        "cc-pVDZ-F12-MP2fit",
        "cc-pVTZ-F12-MP2fit",
        "cc-pVQZ-F12-MP2fit",
        "cc-pCVDZ-F12-MP2fit",
        "cc-pCVTZ-F12-MP2fit",
        "cc-pCVQZ-F12-MP2fit",
        "cc-pVDZ-PP-F12-MP2fit",
        "cc-pVTZ-PP-F12-MP2fit",
        "cc-pVQZ-PP-F12-MP2fit",
        "aug-cc-pwCVDZ-PP/C",
        "aug-cc-pwCVTZ-PP/C",
        "aug-cc-pwCVQZ-PP/C",
        "cc-pwCVDZ-PP/C",
        "cc-pwCVTZ-PP/C",
        "cc-pwCVQZ-PP/C",
    ]
)
df9.columns = ["Auxiliary Correlation Consistent Basis Sets"]
df_qm_basis_sets_I = pd.concat([df1, df2, df3, df4, df5, df6], axis=1)
df_qm_basis_sets_I.fillna("-", inplace=True)
df_qm_basis_sets_I.reset_index(drop=True, inplace=True)
df_qm_basis_sets_II = pd.concat([df7, df8, df9], axis=1)
df_qm_basis_sets_II.fillna("-", inplace=True)
df_qm_basis_sets_II.reset_index(drop=True, inplace=True)
# QM Methods
df1 = pd.DataFrame(
    [
        "MP2",
        "RI-MP2",
        "SCS-MP2",
        "RI-SCS-MP2",
        "OO-RI-MP2",
        "OO-RI-SCS-MP2",
        "MP2-F12",
        "MP2-F12-RI",
        "MP2-F12D-RI",
    ]
)
df1.columns = ["Perturbation Methods"]
df2 = pd.DataFrame(
    [
        "CCSD",
        "CCSD(T)",
        "CCSD-F12",
        "CCSD(T)-F12",
        "CCSD-F12/RI",
        "CCSD-F12D/RI",
        "CCSD(T)-F12/RI",
        "CCSD(T)-F12D/RI",
        "QCISD",
        "QCISD(T)",
        "QCISD-F12",
        "QCISD(T)-F12",
        "QCISD-F12/RI",
        "QCISD(T)-F12/RI",
        "CPF/1",
        "NCPF/1",
        "CEPA/1",
        "NCEPA/1",
        "RI-CEPA/1-F12",
        "MP3",
        "SCS-MP3",
    ]
)
df2.columns = ["Reference Methods"]
df3 = pd.DataFrame(
    [
        "HFS",
        "LDA",
        "VWN",
        "VWN3",
        "PWLDA",
        "BP86",
        "BLYP",
        "OLYP",
        "GLYP",
        "XLYP",
        "PW91",
        "mPWPW",
        "mPWLYP",
        "PBE",
        "RPBE",
        "REVPBE",
        "RPW86PBE",
        "PWP",
    ]
)
df3.columns = ["Local Gradient Corrected Methods"]
df4 = pd.DataFrame(
    [
        "B1LYP",
        "B3LYP",
        "O3LYP",
        "X3LYP",
        "B1P",
        "B3P",
        "B3PW",
        "PW1PW",
        "mPW1PW",
        "mPW1LYP",
        "PBE0",
        "REVPBE0",
        "REVPBE38",
        "BHANDHLYP",
    ]
)
df4.columns = ["Hybrid Methods"]
df5 = pd.DataFrame(
    [
        "TPSS",
        "TPSSh",
        "TPSS0",
        "M06L",
        "M06",
        "M062X",
        "PW6B95",
        "B97M-V",
        "B97M-D3BJ",
        "B97M-D4",
        "PBE0",
        "SCANfunc",
    ]
)
df5.columns = ["Meta GGA / Hybrid Meta GGA Methods"]
df_qm_methods = pd.concat([df1, df2, df3, df4, df5], axis=1)
df_qm_methods.fillna("-", inplace=True)
df_qm_methods.reset_index(drop=True, inplace=True)
# QM/QM2 charge methods
df_qm_qm2_charges = pd.DataFrame(["HIRSHFELD", "CHELPG", "MULLIKEN", "LOEWDIN"])
df_qm_qm2_charges.columns = ["QM / QM2 Charge Methods"]
df_qm_qm2_charges.reset_index(drop=True, inplace=True)

fig, ax = plt.subplots(figsize=(14, 8))
ax.axis("tight")
ax.axis("off")
the_table = ax.table(cellText=df_qm_methods.values, colLabels=df_qm_methods.columns)
the_table.auto_set_font_size(False)
the_table.set_fontsize(15)
the_table.scale(2, 2)
ax.set_title("QM Methods", fontsize=15)
pp = PdfPages("qm_methods.pdf")
pp.savefig(fig, bbox_inches="tight")
pp.close()

fig, ax = plt.subplots(figsize=(14, 8))
ax.axis("tight")
ax.axis("off")
the_table = ax.table(
    cellText=df_qm_basis_sets_I.values, colLabels=df_qm_basis_sets_I.columns
)
the_table.auto_set_font_size(False)
the_table.set_fontsize(15)
the_table.scale(2, 2)
ax.set_title("QM Basis Sets (Part I)", fontsize=15)
pp = PdfPages("qm_basis_sets_i.pdf")
pp.savefig(fig, bbox_inches="tight")
pp.close()

fig, ax = plt.subplots(figsize=(14, 8))
ax.axis("tight")
ax.axis("off")
the_table = ax.table(
    cellText=df_qm_basis_sets_II.values, colLabels=df_qm_basis_sets_II.columns
)
the_table.auto_set_font_size(False)
the_table.set_fontsize(15)
the_table.scale(2, 2)
ax.set_title("QM Basis Sets (Part II)", fontsize=15)
pp = PdfPages("qm_basis_sets_ii.pdf")
pp.savefig(fig, bbox_inches="tight")
pp.close()

fig, ax = plt.subplots(figsize=(14, 8))
ax.axis("tight")
ax.axis("off")
the_table = ax.table(cellText=df_qm2_methods.values, colLabels=df_qm2_methods.columns)
the_table.auto_set_font_size(False)
the_table.set_fontsize(15)
the_table.scale(2, 2)
ax.set_title("QM2 Methods", fontsize=15)
pp = PdfPages("qm2_methods.pdf")
pp.savefig(fig, bbox_inches="tight")
pp.close()

fig, ax = plt.subplots(figsize=(14, 8))
ax.axis("tight")
ax.axis("off")
the_table = ax.table(
    cellText=df_qm_qm2_charges.values, colLabels=df_qm_qm2_charges.columns
)
the_table.auto_set_font_size(False)
the_table.set_fontsize(15)
the_table.scale(2, 2)
ax.set_title("QM/QM2 Charge Methods", fontsize=15)
pp = PdfPages("qm_qm2_charges.pdf")
pp.savefig(fig, bbox_inches="tight")
pp.close()

pdf_list = [
    "qm_methods.pdf",
    "qm_basis_sets_i.pdf",
    "qm_basis_sets_ii.pdf",
    "qm2_methods.pdf",
    "qm_qm2_charges.pdf",
]
merger = PdfMerger()
for pdf in pdf_list:
    merger.append(pdf)
merger.write("orca_methods_basis_sets.pdf")
merger.close()
for i in pdf_list:
    command = "rm -rf " + i
    os.system(command)
