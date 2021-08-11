# MSA

## Summary

For the multiple sequence alignment the following protein sequences from Uniprot were used: Q14643, Q14573, A0A0R4IH08 (Trembl entry, but the Swissprot F1R1L5 is a fragment), F8W4B1, E7FB06, A0A4W3JHJ0, A0A4W3JQ24, A0A4W3GHU7, P29993) and NCBI (XP_032819189.1 (absent in Uniprot). Additionally three NCBI sequences were used: XP_022786205.1 (full-length unlike Trembl A0A2B4T0S0), XP_021323455.1 1 (full-length unlike Trembl F1QW18). These represent the genes from the following species: Homo sapiens (human, crown group), Danio rerio (zebra fish, a boney fish), Callorhinchus milii (Ghost shark, a jawed cartilaginous fish), Petromyzon marinus (sea lamprey, a jawless fish, i.e. diverged before the tetrapseudoploidy event of the jawed fish/tetrapods clade), Drosophila melanogaster (fruit-fly, a protostome invertebrate) and Stylophora pistillata (smooth cauliflower coral, a coral, i.e. a basal animal). Note, a further duplication of ITPR1 occurred in zebrafish.
## Plan

Make a multiple sequence alignment and make a picture in Plotly.

## Methods

Define Uniprot ids:

```python
uniprot_ids = dict(itpr1_human = 'Q14643',
                    itpr2_human = 'Q14571',
                    itpr3_human = 'Q14573',
                    itpr3_zebrafish = 'A0A0R4IH08', #'F1R1L5',
                   # itpr2_zebrafish = 'F1QW18',
                    itpr1a_zebrafish = 'F8W4B1',
                    itpr1b_zebrafish = 'E7FB06',
                    itpr1_shark = 'A0A4W3JHJ0',
                    itpr2_shark = 'A0A4W3JQ24',
                    itpr3_shark = 'A0A4W3GHU7',
                    itpr_coral = 'A0A2B4T0S0',
                    itpr_fly = 'P29993')
```

Define getters:

```python
import requests

def get_fasta(url:str) -> str:
    reply = requests.get(url)
    assert reply.status_code == 200, f'{reply.status_code}: {url}'
    return reply.text

def split_fasta(fasta:str):
    lines = fasta.split('\n')
    header = lines[0].strip()
    data = ''.join([l.strip() for l in lines[1:]])
    return header, data

def get_seq(uniprot_id:str):
    _, seq = split_fasta(get_fasta(f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta'))
    return seq
```

Get seqs and write them:

```python
seqs = [get_fasta(f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta') for uniprot_id in uniprot_ids.values()]
with open('raw.fasta', 'w') as fh:
    fh.writelines(seqs)
```

Check:
```python
{name: len(get_seq(uniprot_id)) for name, uniprot_id in uniprot_ids.items()}
```

Spiken in with three extra sequences.

Align them with muscle --> 'al.fasta'

Correct the names for the fasta.

```python
inverted_uniprot_ids = dict(zip(uniprot_ids.values(), uniprot_ids.keys()))
inverted_uniprot_ids['XP_032819189.1'] = 'itpr_lamprey'
inverted_uniprot_ids['XP_022786205.1'] = 'itpr_coral'
inverted_uniprot_ids['XP_021323455.1'] = 'itpr2_zebrafish'

lines = []
with open('al.fasta', 'r') as fh:
    for line in fh:
        if '>' in line:
            # >sp|P29993|ITPR_DROME
            parts = line.split('|')
            
            lines.append(">"+inverted_uniprot_ids[parts[1]]+'\n')
        else:
            lines.append(line)
with open('al.fixed.fasta', 'w') as fh:
    fh.writelines(lines)
```

And make a dictionary of sequences:

```python
inverted_uniprot_ids = dict(zip(uniprot_ids.values(), uniprot_ids.keys()))
al_seqs = {}
current = 'error'
with open('al.fasta', 'r') as fh:
    for line in fh:
        if '>' in line:
            # >sp|P29993|ITPR_DROME
            parts = line.split('|')
            current = inverted_uniprot_ids[parts[1]]
            al_seqs[current] = ''
        else:
            al_seqs[current] += line.strip()
```

Define the MSA based sequence converters

```python
# 
def get_mapping(seq):
    return [i for i, p in enumerate(seq) if p != '-']

def check_mapping(m, gene_name, t):
    assert ungap(al_seqs[gene_name])[t] == al_seqs['itpr_fly'][m[t]]

def ungap(seq):
    return seq.replace('-','')
    
def convert(query_seq, target_seq, position):
    query_mapping = get_mapping(query_seq)
    target_mapping = get_mapping(target_seq)
    mp = query_mapping[position - 1]
    return target_mapping.index(mp) + 1

def check_covert(query_seq, position):
    convert(query_seq, query_seq, position)
```

Make graphs

```python
import plotly.graph_objects as go
import numpy as np
from collections import Counter


n = len(one)
names = ['itpr3_human', 'itpr3_danio', 'itpr3_shark', 'itpr2_human', 'itpr2_danio', 'itpr2_shark', 'itpr1_human', 'itpr1a_danio', 'itpr1b_danio', 'itpr1_shark','itpr_fly', 'itpr_coral']
def make_fig(al_seqs: Dict[str, str], position: int, names: Optional[List[str]] = None) -> go.Figure:
    """
    This make a fake scatter plot with letters for each position in Plotly.

    * al_seq name ->
    * names, the order wanted for these

    :param al_seqs: aligned seq dictionary name -> seq w/ gaps
    :param position: position to show.
    :param names: order or the seqs to show top to bottom
    :return:
    """
    # --------------------------------------------------
    # prep
    if names is None:
        names = list(al_seqs.keys())
    first = al_seqs[names[0]]
    n = len(first)
    mapping = get_mapping(first)
    position = position - 1
    start = mapping[position] - 10
    mid = mapping[position]
    stop = mapping[position] + 10
    if start < 0:
        forepadding  = ['>'] * abs(start)
        mid += abs(start)
        start = 0
    else:
        forepadding = []
    if stop > n:
        aftpadding  = ['<'] * stop - n
        stop = n
    else:
        aftpadding = []
    # --------------------------------------------------
    # Add letters
    scatters = [go.Scatter(x=np.arange(start, stop + 1),
                           # y=[-i/yzoom] * 20,
                           y=[name] * 21,
                           text=forepadding + list(al_seqs[name][start:stop + 1]) + aftpadding,
                           textposition="middle center",
                           # textfont_size=10,
                           textfont=dict(
                               family="monospace",
                               size=18,
                               color="black"
                           ),
                           name=name,
                           mode='text') for i, name in enumerate(names)
                ]
    # --------------------------------------------------
    # Add red box
    fig = go.Figure(data=scatters)
    fig.add_shape(type="rect",
                  x0=mid - 0.5, x1=mid + 0.5,
                  # y0=0.5/yzoom, y1=(0.5 -len(al_seqs))/yzoom,
                  y0=0 - 0.5, y1=len(names) - 1 + 0.5,  # the len is off-by-one
                  line=dict(color="crimson"),
                  )
    fig.update_shapes(dict(xref='x', yref='y'))
    # --------------------------------------------------
    # Add turquoise squares
    for i in range(start, stop + 1):
        residues = [al_seqs[name][i] for name in names]
        mc = Counter(residues).most_common()[0][0]
        for j, r in enumerate(residues):
            if mc == '-':
                continue
            elif r == mc:
                fig.add_shape(type="rect",
                              x0=i - 0.5, x1=i + 0.5,
                              # y0=(-j -.5)/yzoom, y1=(-j +.5)/yzoom,
                              # y0=names[j], y1=names[j],
                              y0=j - 0.5, y1=j + 0.5,
                              fillcolor="aquamarine", opacity=0.5,
                              layer="below", line_width=0,
                              )
    # --------------------------------------------------
    # Correct
    fig.update_layout(template='none',
                      showlegend=False,
                      title=f'Residue {position + 1}',
                      xaxis={
                          'showgrid': False,
                          'zeroline': False,
                          'visible': False,
                          # 'range': [start,stop],
                      },
                      yaxis={
                          'showgrid': False,
                          'zeroline': False,
                          'visible': True,
                          'autorange': 'reversed'
                      }
                      )
    # --------------------------------------------------
    return fig
```

Make for each position:

```python
fig = make_fig(2473)
#fig.write_image(f"MSA/fig_resi{position}.svg")
fig.show()
```

## Equivalent residues

Find what the positions are of the ITPR1 variants.
There is the detail that NP_001161744.1 (2743 aa) is isoform 3. So it needs converting too.

```python

from Bio import pairwise2
alignment = pairwise2.align.globalxx(isoform3, ungap(al_seqs['itpr1_human']))[0]

shitty_one = alignment[0]
correct_one = alignment[1]

import re
for m in "p.Gly2539Arg", "p.Lys2596del", "p.Glu2094Gly":
    position = int(re.search(r'(\d+)', m).group(1))
    one = al_seqs['itpr1_human']
    three = al_seqs['itpr3_human']
    corposition = convert(shitty_one, correct_one, position)
    transposition = convert(one, three, corposition)
    print('mutation', m)
    print(f'Stated {position}, corrected {corposition}, to 3 {transposition}')
    print(f'In ITPR1...', ungap(one)[corposition - 1], ungap(one)[corposition - 5: corposition + 5])
    print(f'In ITPR3...', ungap(three)[transposition - 1], ungap(three)[transposition - 5: transposition + 5])
```

    mutation p.Gly2539Arg
    Stated 2539, corrected 2554, to 3 2473
    In ITPR1... G GLRSGGGVGD
    In ITPR3... G GLRNGGGVGD
    mutation p.Lys2596del
    Stated 2596, corrected 2611, to 3 2530
    In ITPR1... K EKQKKEEILK
    In ITPR3... K EKQKKEEILK
    mutation p.Glu2094Gly
    Stated 2094, corrected 2109, to 3 2005
    In ITPR1... E LAIMESRHDS
    In ITPR3... E LALMESRHDS
    