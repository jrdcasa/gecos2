# Tutorial 1: Prepare conformers from a PDB file.
----

**Folder**        : 01-EVOH_model01_aoki/04-CONF_OPENBABEL_DRAGO

**Input file**    : ../01-CONF_RDKIT_AOKI/Model01_EVOH_ini.pdb

**Python script** : ../004-Model01-EVOH.py

**Json script**   : ../004-Model01-EVOH.json

**Video** : TODO

This tutorial shows you how to create conformers from a PDB file. The PDB files must contain the section "CONECT". 

Start the python environment in which GeCos is installed and execute the command **gecos_gui** (see [installation section](./02-installation.md)). The gui allows one to create a python script to run GeCos. 

## Conformers with RDKIT algotithm

The steps to generate conformers with RdKit are approximately the following:

    * Generate the number of conformers set in the user interface.
    * Minimise the conformers with a force field 
    * Align the minimised conformers and group the conformers based on the RMSD of the heavy atoms. 
    * The conformers resulting from the last step are those that are optimised in the QM program.