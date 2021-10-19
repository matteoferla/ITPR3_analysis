# ITPR3_analysis
In silico analyses of ITPR3 variants.
This is the support data for a [Michelanglo page](https://michelanglo.sgc.ox.ac.uk/) —page yet to be released.

## Models

Inositol 1,4,5-trisphosphate receptor type 3 ([Uniprot: Q14573](https://www.uniprot.org/uniprot/Q14573))
is a calcium channel activated by inositol 1,4,5-trisphosphate.

In [Paknejad & Hite, 2018](https://www.nature.com/articles/s41594-018-0089-6), a variety of structures were solved
by cryoEM in various states. 

|    | name                 | emd      | pdb   |  Zn   |  I3P  |  Ca   | Presumed state |
|---:|:---------------------|:---------|:------|:------|:------|:------|:---------------------|
|  0 | hIP3R3 apo           | EMD-7978 | 6DQJ  | 1     |  0    | 0     |Closed                |
|  1 | hIP3R3 IP3 class 1   | EMD-7981 | 6DQN  | 1     |  1    | 0     |Pre-activation 1:1:1:1|
|  2 | hIP3R3 IP3 class 2   | EMD-7984 | 6DQV  | 1     |  1    | 0     |Pre-activation 2:2:2:2| 
|  3 | hIP3R3 IP3 class 3   | EMD-7983 | 6DQS  | 1     |  1    | 0     |Pre-activation 2:2:2:1|
|  4 | hIP3R3 IP3 class 4   | EMD-7986 | 6DQZ  | 1     |  1    | 0     |Pre-activation 2:2:1:1|
|  5 | hIP3R3 IP3 class 5   | EMD-7987 | 6DR0  | 1     |  1    | 0     |Pre-activation 2:1:1:1|
|  6 | hIP3R3 Ca2+ bound    | EMD-7988 | 6DR2  | 1     |  0    | 2     |Inhibited             |
|  7 | hIP3R3 low IP3–Ca2+  | EMD-7991 | 6DRA  | 1     |  0    | 2     |Inhibited             | 
|  8 | hIP3R3 high IP3–Ca2+ | EMD-7994 | 6DRC  | 1     |  1    | 2     |Inhibited             |

(The class 1 and class 2 have the ARM2 domain in different positions, class 1 is similar to Apo, while 2 isn't,
all other classes are mixes of class 1 and class 2 states).

Do note that the models exceed 99,999 atoms making the PDB format problematic.
As a consequence the minimised models are present in the [folder structures](structures) as
both `pdb` and `mmtf`.

| definition    |   from |   to |
|:--------------|-------:|-----:|
| cytoplasmic   |      1 | 2202 |
| transmembrane |   2206 | 2671 |
| Ins145_P3_rec |      3 |  230 |
| MIR           |    233 |  433 |
| RYDR_ITPR     |    473 |  670 |
| RYDR_ITPR     |   1175 | 1334 |
| RIH_assoc     |   1864 | 1974 |
| Ion_trans     |   2206 | 2528 |

## Calculations

Pyrosetta was used.
Base script are in [my pyrosetta_scripts repo](https://github.com/matteoferla/pyrosetta_scripts).
Code used is in [code notes](code.md).

Each of the structures was minimised using LocalRelax mover and against the respective electron density.

In light of the absence of key membrane facing residues, the membrane framework was not used.

Given a mutation, the 12 Å neighbourhood of the wild type was energy minimised (resulting in the ∆G of the wild type),
the mutation introduced and minimised again.

A final 5 cycles of dihedral Relax per chain was initially considered by not done (5 CPU days and insufficient)
hence the 12 Å neighbourhood thorough minimisation was done for each variant.

The ref2015 scorefunction was used with a 2.94 scaling to convert REU to kcal/mol. 

![heatmap](figures/heatmap3-01.png)

### Homozygous

![quad](figures/quad-01.png)

Three residues are close enough to affect the scores of that same residue on other chains.

### gnomAD

Residues in regions of interest (see [code notes](code.md)):

p.Ala18Ser, p.Asn23Ser, p.Asn23Lys, p.Ile184Ser, p.Trp219Ter, p.Arg270His, p.Asp742Glu, p.Arg743Cys, p.Val1121Met, p.Thr1945Ile, p.Glu1946ValfsTer8, p.Arg2471His, p.Gly2473Ser, p.Gly2477Ala, p.Asp2478Asn, p.Leu2511Val, p.Val2515Leu, p.Asp2582His, p.Asp2582Asn, p.Tyr2583Asn

Top 10 most destabilizing, with Al's cutoff of 5e-5.


### lipid
A side investigation was done to see if the pocket formed by 2527, 2524 and 2357 could bind an organophosphate
(e.g. phospholipid). The pocket can indeed favourably accomodate an organophosphate (ethylphosphate as a proxy),
with a ∆∆G_binding of –17 kcal/mol, resulting from the following contributions:


However, in the density map it appears suitable for a water. While in [Baker et al. 2021](https://www.nature.com/articles/s42003-021-02156-4) the cryoEM structure of ITPR1 with
lauryl maltose neopentyl glycol (LMNG, residue name: PLX) was solved (PDB:7LHF), in this model the lipids are close
but there is not a phosphate in that exact pocket, so this investigation was halted.

