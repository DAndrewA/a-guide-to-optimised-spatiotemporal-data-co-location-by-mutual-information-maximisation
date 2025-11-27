"""Author: Andrew Martin
Creation date: 27/11/25

Script handling functionality to reporoduce Table 2 showing results from copula analysis.
"""

import sys
sys.path.insert(1, "../")
from common.handle_MI_datasets import (
    K
)
from common.handle_vcfs import vcfs_per_parametrisation, generate_confusion_matrix, generate_masks, PARAMETRISATION_print_names
from common.copula import BivariateCopula


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ecdf


from dataclasses import dataclass
from typing import Self, Any
import numpy as np
import xarray as xr


@dataclass(kw_only=True)
class TableEntry:
    value: Any
    bold: bool

    def __str__(self):
        return f"{self.value}"

@dataclass(kw_only=True)
class NumericTableEntry(TableEntry):
    sf: int

    def __str__(self):
        return "".join([
            "$",
            f"{r'\mathbf{' if self.bold else ''}",
            f"{self.value:{f'.{self.sf}f'}}",
            f"{r'}' if self.bold else ''}",
            "$"
        ])

@dataclass(kw_only=True)
class TableRow:
    plabel: TableEntry
    R_km: NumericTableEntry
    tau_s: NumericTableEntry
    accuracy: NumericTableEntry
    cmin: NumericTableEntry
    cmax: NumericTableEntry
    c11: NumericTableEntry
    RMSD: NumericTableEntry

    def __str__(self):
        return " & ".join([
            str(field)
            for field in (
                self.plabel, 
                self.R_km, 
                self.tau_s, 
                self.accuracy, 
                self.cmin, 
                self.cmax, 
                self.c11, 
                self.RMSD
            )
        ])

    @classmethod
    def from_values(cls, plabel, R_km, tau_s, accuracy, cmin, cmax, c11, RMSD) -> Self:
        return cls(
            plabel = TableEntry(value=plabel, bold=False),
            R_km = NumericTableEntry(value=R_km, bold=False, sf=1),
            tau_s = TableEntry(value=tau_s, bold=False),
            accuracy = NumericTableEntry(value=accuracy, bold=False, sf=3),
            cmin = NumericTableEntry(value=cmin, bold=False, sf=2),
            cmax = NumericTableEntry(value=cmax, bold=False, sf=2),
            c11 = NumericTableEntry(value=c11, bold=False, sf=2),
            RMSD = NumericTableEntry(value=RMSD, bold=False, sf=2),
        )

CAPTION = "\n    ".join(
r"""
Results from the computation of copulae comparing VCF values between ATL09 and Cloudnet data for all the tested parametrisations.
The accuracy of the agreement between the VCF retrievals for the categories nc, pc and tc (see Figure \ref{fig::copulae}) is given as $\text{ACC}$.
$c_\text{min}$ is the minimum copula density for the given parametrisation, $c_\text{max}$ is the maximum achieved copula density, and $c(1,1)$ is the discretised tail dependence of the copula.
$\text{RMSD}$ is the root mean squared difference of the copula density from the independence copula density.
Values in bold indicate the best parametrisation for the given metric (the notion of best being defined in the text).
""".split("\n")
)

table_pre = "\n".join([
    r"""
%t
\begin{table*}[t]
\caption{""",
    CAPTION,
    r"""}
\label{table::copulae}
\begin{tabular}{lccccccr}
\tophline
    parametrisation & $R$ (\unit{km}) & $\tau$ (\unit{s}) & $\text{ACC}$ & $c_\text{min}$ & $c_\text{max}$ & c(1,1) & $\text{RMSD}$ \\
\middlehline
"""
])

table_post = r"""
\bottomhline
\end{tabular}
\belowtable{} % Table Footnotes
\end{table*}
"""


# load the VCF data
vcfs_per_p = vcfs_per_parametrisation()

# remove unwanted parametrisations
del vcfs_per_p["P_01"], vcfs_per_p["P_10"]

table_rows = list()

for plabel, vcfs in vcfs_per_p.items():
    # compute the accuracy of the confusion matrix
    confusion_matrix = generate_confusion_matrix(vcfs)
    normalised_confusion_matrix = confusion_matrix / np.sum(confusion_matrix)
    accuracy = np.diag(normalised_confusion_matrix).sum()

    # generate the masks that define the (pc,pc) subset of interest 
    (_, pc_atl09, _) = generate_masks(
        data = (data_atl09:=vcfs.vcf_atl09.data)
    )
    (_, pc_cloudnet, _) = generate_masks(
        data = (data_cloudnet:=vcfs.vcf_cloudnet.data)
    )
    mask_non_degenerate = pc_atl09 & pc_cloudnet

    copula = BivariateCopula.generate(
        data_X = data_atl09[mask_non_degenerate].flatten(),
        data_Y = data_cloudnet[mask_non_degenerate].flatten(),
    )

    table_rows.append(
        TableRow.from_values(
            plabel = PARAMETRISATION_print_names[plabel],
            R_km = 0.,
            tau_s = 0,
            accuracy = accuracy,
            cmin = copula.cmin,
            cmax = copula.cmax,
            c11 = copula.c11,
            RMSD = copula.RMSD,
        )
    )


# identify the best value across the table rows for accuracy,...,RMSD
best_accuracy = table_rows[0]
best_cmin = table_rows[0]
best_cmax = table_rows[0]
best_c11 = table_rows[0]
best_RMSD = table_rows[0]
for row in table_rows[1:]:
    if row.accuracy.value > best_accuracy.accuracy.value: best_accuracy = row
    if row.cmin.value < best_cmin.cmin.value: best_cmin = row
    if row.cmax.value > best_cmax.cmax.value: best_cmax = row
    if row.c11.value > best_c11.c11.value: best_c11 = row
    if row.RMSD.value > best_RMSD.RMSD.value: best_RMSD = row

best_accuracy.accuracy.bold = True
best_cmin.cmin.bold = True
best_cmax.cmax.bold = True
best_c11.c11.bold = True
best_RMSD.RMSD.bold = True

table_string = "\n".join([
    table_pre,
    "\n".join([
        str(row)
        for row in table_rows
    ]),
    table_post,
])

with open("table-2.txt", "w") as f:
    f.write(table_string)

