# ITPR3_analysis
In silico analyses of ITPR3 variants

## Models

Inositol 1,4,5-trisphosphate receptor type 3 ([Uniprot: Q14573](https://www.uniprot.org/uniprot/Q14573))
is a calcium channel activated by inositol 1,4,5-trisphosphate.

In [Paknejad & Hite, 2018](https://www.nature.com/articles/s41594-018-0089-6), a variety of structures were solved
by cryoEM in various states. 

|    | name                 | emd      | pdb   |
|---:|:---------------------|:---------|:------|
|  0 | hIP3R3 apo           | EMD-7978 | 6DQJ  |
|  1 | hIP3R3 IP3 class 1   | EMD-7981 | 6DQN  |
|  2 | hIP3R3 IP3 class 2   | EMD-7984 | 6DQV  |
|  3 | hIP3R3 IP3 class 3   | EMD-7983 | 6DQS  |
|  4 | hIP3R3 IP3 class 4   | EMD-7986 | 6DQZ  |
|  5 | hIP3R3 IP3 class 5   | EMD-7987 | 6DR0  |
|  6 | hIP3R3 Ca2+ bound    | EMD-7988 | 6DR2  |
|  7 | hIP3R3 low IP3–Ca2+  | EMD-7991 | 6DRA  |
|  8 | hIP3R3 high IP3–Ca2+ | EMD-7994 | 6DRC  |

These were used for scoring.

## Calculations

Pyrosetta was used.
Base script are in [my pyrosetta_scripts repo](https://github.com/matteoferla/pyrosetta_scripts).
Code used is in [code notes](code.md).

Namely, the structures were minimised using local relax against the density map.
A final 5 cycles of dihedral Relax per chain was initially considered by not done (5 CPU days)
instead a 12 Å neighbourhood thorough minimisation was done for each variant.

