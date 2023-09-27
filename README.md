# Proton Disorder

## Background
Ice structures obtained via X-ray Crystallography studies can only parse out the position of the oxygen atoms in the structure. Therefore, it is necessary to come up with a method for finding the position of the hydrogen atoms in the crystal structure, such that you have a water crystal.

The position of the hydrogens can be found using some basic rules:
1) Each oxygen atom will have 4 neighboring oxygen atoms. The 2 hydrogens for each oygen will point at another oxygen.
2) If you have a 'link' between two oxygens, only one hydrogen can reside in that link.
3) The dipole moment for the entire structure should be ~0.

To accomplish this, `proton-disorder` will parse a x-ray crystallography structure with oxygen atom positions, find the 'links' between each oxygen pair, and then place a hydrogen in each one of these links. It will then make sure that one of the oxygens 'owns' this hydrogen, such that each oxygen has two hydrogens associated with it.

`proton-disorder` will then follow a Monte-Carlo algorithm to move the position of the hydrogens around an oxygen, adjust the hydrogens again so that every oxygen has two hydrogens, and then calculate the dipole moment of the ice. If the dipole moment is lower than the previous iteration, it will keep that structure and keep shaking the hydrogen bonding until the dipole moment is beneath the threshold.

## Installation

First, create a virtual python environment using `conda` or `venv`

### Conda
``` bash
$ conda create -n proton-disorder
$ conda activate proton-disorder
```

### Venv
_Instructions have not been confirmed_
``` bash
$ python3 -m venv proton-disorder
$ source /path/to/venv/bin/activate
```

Then, clone this repository:
``` bash
$ git clone https://github.com/BoiseState-MATERIALab/proton-disorder.git
```

Navigate into the cloned directory, and run the following command to install `proton-disorder` and dependencies into your virtual environment:
``` bash
$ python3 -m pip install .
```

## Usage

To view help message, use the following command:
```bash
$ python3 -m proton_disorder -h
```

Output
```
usage: Proton Disorder [-h] [-i INPUT_FILE] [-o OUTPUT_FILE] [-d DIPOLE_TARGET] [-n NSHAKES] [-t]

Generates a proton disordered ice structure

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file INPUT_FILE
                        Name of input xyz file
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        Name for output file
  -d DIPOLE_TARGET, --dipole-target DIPOLE_TARGET
                        Target dipole for hydrate
  -n NSHAKES, --nshakes NSHAKES
                        Number of shakes to do per round
  -t, --tip3p           Use TIP3P (Not recommended)

nerd
```

The only requred parameter is the input filename (`-i` or `--input-file`). The input file should be formatted in `csv` or `xyz` format. The only difference between the two are delimination (`csv` uses commas, `xyz` uses tabs). An example of each can be found in the [input](input) directory. 

Currently the resulting ice structure can be outputted as an `xyz`, `csv`, or `pdb` file. Just use the `-o` or `--output-file` to specifiy the output file path. Specify the file format in the output filename.

The recommended dipole target is 0.1 or lower. It is recommended to also use a value for nshakes that is ~1/10 the number of atoms in your system, which is calculated for you be default.

TIP3P water can also be written for the ice output, however this is not recommended. For most applications, TIP4P would be preferred.
