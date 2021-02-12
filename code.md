## Downloads

IDs from Paknejad & Hite, 2018

    raw = 'hIP3R3 apo (EMD-7978, PDB 6DQJ)	hIP3R3 IP3 class 1 (EMD-7981, PDB 6DQN)	hIP3R3 IP3 class 2 (EMD-7984, PDB 6DQV)	hIP3R3 IP3 class 3 (EMD-7983, PDB 6DQS)	hIP3R3 IP3 class 4 (EMD-7986, PDB 6DQZ)	hIP3R3 IP3 class 5 (EMD-7987, PDB 6DR0)	hIP3R3 Ca2+ bound (EMD-7988, PDB 6DR2)	hIP3R3 low IP3–Ca2+ (EMD-7991, PDB 6DRA)	hIP3R3 high IP3–Ca2+ (EMD-7994, PDB 6DRC)'
    
    import re
    import pandas as pd
    
    models = [dict(name = re.match('(.*)\(', entry).group(1),
                   emd = re.search('(EMD-\d+)', entry).group(1),
                   pdb = re.search('PDB (\w{4})', entry).group(1),
                  ) for entry in raw.split('\t')]
    model_table = pd.DataFrame(models)

Data retrieval: EM map

    import shutil
    import urllib.request as request
    from contextlib import closing
    
    def download_map(code: str):
        ftp_path = f'ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/{code}/map/emd_{code.replace("EMD-","")}.map.gz'
        file_path = f'{code}.map.gz'
        with closing(request.urlopen(ftp_path)) as r, open(file_path, 'wb') as f:
            shutil.copyfileobj(r, f)
            
Let's download cifs this time.
            
    import requests

    def download_pdb(code: str):
        """
        Download CIF
        """
        http_path = f'http://www.ebi.ac.uk/pdbe/static/entry/download/{code.lower()}-assembly-1.cif.gz'
        file_path = f'{code}.cif.gz'
        r = requests.get(http_path, stream=True)
        if r.status_code == 200:
            with open(file_path, 'wb') as f:
                r.raw.decode_content = True
                shutil.copyfileobj(r.raw, f)
            return True
        else:
            return False
            
OMP aligned

    import requests
    
    def download_opm(code: str):
        """
        Download OPM PDB
        """
        http_path = f'https://opm-assets.storage.googleapis.com/pdb/{code.lower()}.pdb'
        file_path = f'{code}_OMP.pdb'
        r = requests.get(http_path, stream=True)
        if r.status_code == 200:
            with open(file_path, 'wb') as f:
                r.raw.decode_content = True
                shutil.copyfileobj(r.raw, f)
            return True
        else:
            return False
            
## Any ligands?

    import pymol2
    
    for i, row in model_table.iterrows():
        pdb_filename = row.pdb+'.cif'
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(pdb_filename, 'model')
            ligands = {atom.resn for atom in pymol.cmd.get_model('not polymer').atom}
            print(ligands)
            
Gives 'ZN', 'I3P' and 'CA'.

Smiles from PDB needs correcting as I3P is deffo not neutral

    from rdkit import Chem
    #Chem.MolFromSmiles('O[C@@H]1[C@H](O)[C@@H](O[P](O)(O)=O)[C@H](O[P](O)(O)=O)[C@@H](O)[C@@H]1O[P](O)(O)=O')
    smiles = 'O[C@@H]1[C@H](O)[C@@H](O[P]([O-])([O-])=O)[C@H](O[P]([O-])([O-])=O)[C@@H](O)[C@@H]1O[P]([O-])([O-])=O'
    mol = Chem.MolFromSmiles(smiles)
    
![I3P_deproton](I3P_deproton.png)

    from rdkit_to_params import Params
    
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(pdb_filename, 'model')
        pymol.cmd.remove('not resn I3P')
        resi = {atom.resi for atom in pymol.cmd.get_model('not polymer').atom}
        assert len(resi) == 1
        pymol.cmd.save('I3P.pdb')
        
    params = Params.from_smiles_w_pdbfile(pdb_file='I3P.pdb',
                                            smiles=smiles,
                                            name='I3P',
                                            proximityBonding=True)
    params.dump('I3P.params')

Checking it is okay:

    import nglview

    view = nglview.show_rosetta(params.test())
    view.add_representation('hyperball', '*')
    view 
     
## Minimise

Regular boilerplate 

    import pyrosetta
    from init_helper import make_option_string, get_logger
    
    # capture to log
    logger = get_logger()
    
    extra_options= make_option_string(no_optH=False,
                                      ex1=None,
                                      ex2=None,
                                      mute='all',
                                      ignore_unrecognized_res=False, # raise error!
                                      load_PDB_components=False,
                                      ignore_waters=False)
    
    pyrosetta.init(extra_options=extra_options)
    
Some functions:
    
    from typing import *

    def get_local_scorefxn() -> pyrosetta.ScoreFunction:
        ## local scorefxn w/ ED
        # ---------- Local weights ------------------
        # these are mostly the same except for the line:
        # <Set scale_sc_dens_byres="E:0.56,D:0.56,R:0.76,K:0.76,M:0.76,C:0.81,Q:0.81,H:0.81,N:0.81,T:0.81,S:0.81,Y:0.88,W:0.88,A:0.88,F:0.88,P:0.88,I:0.88,L:0.88,V:0.88"/>
        # which is an utter mystery.
        weights = {"cart_bonded_length": 0.5,
                   "cart_bonded_torsion": 0.5,
                   "cart_bonded_angle": 1.0,
                   "pro_close": 0.0,
                   "fa_sol":0.0,
                   "elec_dens_fast": 30, # <-- ED
                   "rama": 0.0,
                   "rama_prepro": 1.0}
        scorefxn_local = pyrosetta.create_score_function('ref2015_cart')
        stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
        for name, value in weights.items():
            scorefxn_local.set_weight(stm.score_type_from_name(name), value)
        return scorefxn_local
    
    def get_pose(pdb_filename: str,
                 params_filenames: Optional[List[str]]) -> pyrosetta.Pose:
        pose = pyrosetta.Pose()
        if params_filenames and isinstance(params_filenames, pyrosetta.rosetta.utility.vector1_string):
            pyrosetta.generate_nonstandard_residue_set(pose, params_filenames)
        if params_filenames and isinstance(params_filenames, list):
            params_filenames2 = pyrosetta.rosetta.utility.vector1_string()
            params_filenames2.extend(params_filenames)
            pyrosetta.generate_nonstandard_residue_set(pose, params_filenames2)
        else:
            pass
        pyrosetta.rosetta.core.import_pose.pose_from_file(pose, pdb_filename)
        return pose
        
    def prep_ED(pose: pyrosetta.Pose, map_filename: str) -> pyrosetta.rosetta.core.scoring.electron_density.ElectronDensity:
        # rmsd & ED fit
        rmsd = pyrosetta.rosetta.core.simple_metrics.metrics.RMSDMetric(pose)
        ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap(map_filename)
        initial_fit = ED.matchPose(pose)
        # This is redundant with ED.
        map_mover = pyrosetta.rosetta.protocols.cryst.LoadDensityMapMover(map_filename)
        map_mover.apply(pose)
        # This is redundant with map
        sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
        sdsm.apply(pose)
        return ED
    
    def local_relax(pose: pyrosetta.Pose, 
                    scorefxn_local: Optional[pyrosetta.ScoreFunction]=None) -> None:
        if scorefxn is None:
            scorefxn = get_local_scorefxn()
        #<LocalRelax name="local_rlx" scorefxn="dens" max_iter="100" ncyc="1" ramp_cart="0" K="16" nexp="2"/>
        relax = pyrosetta.rosetta.protocols.relax.LocalRelax()
        relax.set_sfxn(scorefxn_local)
        relax.set_K(16)
        relax.set_max_iter(100)
        relax.set_ncyc(3)
        relax.set_nexp(2)
        relax.apply(pose)
    
    def chainwise_relax(pose: pyrosetta.Pose, 
                        scorefxn: Optional[pyrosetta.ScoreFunction]=None, 
                        cycles:int=5) -> None:
        if scorefxn is None:
            scorefxn = pyrosetta.get_fa_scorefxn()
        for chain_i in range(1, pose.num_chains()+1):
            chain_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain_i)
            chain_vector = chain_sele.apply(pose)
            movemap = pyrosetta.MoveMap()
            movemap.set_bb(allow_bb=chain_vector)
            movemap.set_chi(allow_chi=chain_vector)
            relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
            relax.set_movemap(movemap)
            relax.apply(pose)
    
Minimise:

    for i, row in model_table.iloc[::-1].iterrows():
        header(row.pdb)
        # ------- files -------------- 
        pdb_filename = row.pdb+'.cif'
        local_filename = f'{row.pdb}.local.pdb'
        map_filename = row.emd+'.map'
        # ------- scorefxns -------------- 
        scorefxn_local = get_local_scorefxn()
        scorefxn = pyrosetta.get_fa_scorefxn()
        elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
        scorefxn.set_weight(elec_dens_fast, 30)
        # --------- Load ------------------
        if os.path.exists(local_filename):
            print('Local already done.')
            pose = get_pose(local_filename, params_filenames=['I3P.params'])
        else:
            pose = get_pose(pdb_filename, params_filenames=['I3P.params'])
        ED = prep_ED(pose, map_filename)
        # --------- Local ------------------
        if not os.path.exists(local_filename):
            local_relax(pose, scorefxn_local)
            print(scorefxn_local(pose), 
                  scorefxn(pose), 
                  rmsd.calculate(pose), 
                  ED.matchPose(pose), 
                  flush=True)
            pose.dump_scored_pdb(local_filename, scorefxn_local)
        # ------- Regular ------------------
        # relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
        # relax.apply(pose)
        # ------- Per Chain ------------------
        chainwise_relax(pose, scorefxn, 5)
        pose.dump_scored_pdb(f'{pdb_name}.lrc.pdb', scorefxn)
        print(scorefxn_local(pose), 
                  scorefxn(pose), 
                  rmsd.calculate(pose), 
                  ED.matchPose(pose), 
                  flush=True)
                  
Due to the problem of size, membrane was not activated.

## Scoring

    from score_mutants import MutantScorer

    scoresx = []
    
    for pdb in model_table.pdb:
        model = MutantScorer.from_file(modelname=pdb,
                                       filename=f'{pdb}.local.pdb',
                                       params_filenames=['I3P.params'])
        model.scorefxn = pyrosetta.create_score_function('ref2015')
        model.strict_about_starting_residue = True
        data = model.score_mutations(variants,
                                    chain='A',
                                    interfaces=(),
                                    preminimise=True,
                                    distance=12,
                                    cycles=5)
        scoresx.append(pd.DataFrame(data))
    
    scores = pd.concat(scoresx)

    scores.to_csv('scores.csv')
    
## Plot

long to wide:

    wide = pd.DataFrame(scores, columns=['model','mutation','complex_ddG'])\
                      .pivot(index='model', columns='mutation', values='complex_ddG')\
                      .round(1)
                      
heatmap
                      
    import plotly.graph_objects as go

    names = dict(zip(model_table.pdb, model_table.name.apply(lambda x: x.strip())))
    
    fig = go.Figure(data=go.Heatmap(
                                    y=wide.columns,
                                    x=wide.index.to_series().apply(lambda x: names[x]),
                                    z=wide.values.transpose(),
                                    colorscale = [(0,"white"),(0.2,"white"), (1,"red")],
                                    zmin=0, zmax=10,
                                    colorbar = {'title': {'text': '∆∆G<br>[kcal/mol]'}}
                                    ),
                   layout=dict(title = {'text': '∆∆G of variants in different structural models'},
                   yaxis = {'title': {'text': 'Variant'}},
                   xaxis = {'title': {'text': 'Model'}},
                   
                ))
    fig.show()
    