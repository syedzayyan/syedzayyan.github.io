+++
title="Docking Preparation Shenanigans in Autodock/Autodock-GPU/Vina"
date="2023-09-19T23:00:01.000Z"

[taxonomies] 
tags = ["tutorials"]
+++

Autodock likes its protein and ligand file in PDBQT formats, which is sort of a pain to do if you are on a MacOS machine. Thus everything from here is Linux related. I have been very lucky to get access to an HPC, but you should be able to do this on a Raspberry Pi or even PythonAnywhere as they have access to bash! 

### Packages

The pip commands to install some of the packages.

```bash
pip install biopython
pip install meeko
```

### Downloading the protein and processing it.

Download the protein. 5N2S can be replaced by your PDB code

```bash
wget http://files.rcsb.org/download/5N2S.pdb
```

In our case the receptor is attached to a cytochrome and I would very much like it gone. I don't really need any of the cytochrome fluff because that is way outside the ligand binding pocket. We hit a problem again. Usually the cytochrome kind of things are in a different chain. Here there are not. You can use the 3D viewer of the RCSB page and figure which residues are worth removing, and then come back to this code.

```python
from Bio import PDB
from Bio.PDB import PDBIO

pdbio = PDBIO()
parser = PDB.PDBParser(QUIET=True)

structure = parser.get_structure("my_structure", "./5N2S.pdb")
list_of_residues = list(range(1002, 1110 + 1)) # Residue IDs worth removing

for chain in structure.get_chains():
    for res in chain.get_residues():
        if res.id[1] in list_of_residues:
            chain.detach_child(res.id)

pdbio.set_structure(structure)
pdbio.save('./5N2S_rem_cyt.pdb')
```

If you double click on any of the residues, the ID comes up on the right.

![RSCB Viewer](/img/blog-img/docking/pdbviewer_rscb.png)

### Remove ligands and stuff and only keep proteins

grep is a good text editing magic thing in unix-based systems.

```bash
grep ATOM 5N2S.pdb > 5N2S_rem.pdb
```



### Download ADFR Suite

We need the ADFR suite to automate a few things that I don't want to break.

```bash
wget -O ADFRsuite_x86_64Linux_1.0.tar.gz https://ccsb.scripps.edu/adfr/download/1038/
```

Sadly enough, Autodock's installation thingy does not work for me on CentOS and it's just easier to download the tarball file and install from there. So that's what we are doing!

After than unzip the file, cd into the uncompressed directory and just run the installation script. At the end of the installation script, there is an alias that you can add to your terminal which makes life easier. All those commands in order :).

```bash
tar zxvf ADFRsuite_x86_64Linux_1.0.tar.gz
cd ADFRsuite_x86_64Linux_1.0
./install.sh
```

Everyone's alias will look different. DON'T Just copy mine. It won't break or anything just you know!

```bash
export PATH=/home/szm25/code mphil/ADFRsuite-1.0/ADFRsuite x86_64Linux 1.0/bin:$PATH
```
I am adding hydrogens because my structure probably does not come with hydrogens. I am learning as well so if you think mine does not drop me and email.

```bash
prepare_receptor -r 5N2S_rem.pdb -o 5N2S.pdbqt -A "hydrogens"
```

### Processing Ligands from SMILES

I use a nifty little Python package called [Meeko](https://github.com/forlilab/Meeko) and here's the full code for that. The function takes in a SMILES string and spits out a PDBQT file in your desired place.

```python
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from rdkit.Chem import AllChem


def make_pdbqt_file_from_smi(smi, file_name, working_dir):
    # Turns the SMILES string into rdkit molecules and protonates it
    try:
        preparator = MoleculePreparation()
        lig = rdkit.Chem.MolFromSmiles(smi)
        protonated_lig = rdkit.Chem.AddHs(lig)
    except:
        print(f"{file_name} SMILES Tucked")

    # Tries to take the rdkit molecule and embed it into 3D space but it's not always for highly flexible compounds.
    # You could try to increase the number of tries but at some point do you want to wait 10 hours for one molecule?
    # The 3D ligand now is turned into PDBQT string by meeko and then churned into a file
    try:
        rdkit.Chem.AllChem.EmbedMolecule(protonated_lig)
        mol_setups = preparator.prepare(protonated_lig)

        for setup in mol_setups:
            pdbqt_string = PDBQTWriterLegacy.write_string(setup)
            with open(f"{working_dir}/{file_name}.pdbqt", "w") as text_file:
                text_file.write(pdbqt_string[0])
    except:
        print(smi, file_name)
```

Now all you need is to run a docking program. There are many tutorials:

- [GNINA](https://colab.research.google.com/drive/1QYo5QLUE80N_G28PlpYs6OKGddhhd931?usp=sharing#scrollTo=WctyMpdMluFN)
- [Autodock GPU](https://www.kaggle.com/code/syedzayyanmasud/adora-docks-duds/settings?scriptVersionId=141909183)
- [Vina](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html#preparing-the-receptor)
- [Smina](https://projects.volkamerlab.org/teachopencadd/talktorials/T015_protein_ligand_docking.html)

I may add one of my own to the pile someday!

Have fun 😊