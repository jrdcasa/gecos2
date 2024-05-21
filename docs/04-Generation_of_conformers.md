# Generation of conformers

Gecos implements two methods for generating conformers:

* **RDKIT**: used the algorithm implemented in [RDKIT toolkit](https://www.rdkit.org/). This toolkit generate conformers using a method called [distance geometry](https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470125823.ch6). The options are:

    * **Number of conformers** _(default: 2000)_: The number of conformers initially generated.

    * **Minimize iterations** _(default 1000)_: Number of itererations for the force field minimization.

    * **Cutoff RMSD Cluster QM(A)** _(default: 0.5)_: RMSD to group conformers after the QM optimizations. 

    * **maxattemps** _(default: 1000)_:  Maximum number of attempts to try embedding
    * **pruneRmsThresh** _(default: -0.01)_: Retain only the conformations out of ‘numConfs’ after embedding that are at least this far apart from each other.  RMSD is computed on the heavy atoms. Pruning is greedy; i.e. the first embedded conformation is retained and from then on only those that are at least pruneRmsThresh away from all retained conformations are kept. The value by default keep this parameter without effect.

    * **cluster threshold**: RMSD criteria to group atoms in a cluster after FF minimization. RMDS is computed on the heavy atoms. Before to calculate the RMSD, the conformers are alligned. 

    * **enforceChirality** : enforce the correct chirality if chiral centers are present. It does not have effect in this version. 

    * **useExpTorsionAnglePrefs** : impose experimental torsion angle preferences. It does not have effect in this version.

    * **useBasicKnowledge** : impose basic knowledge such as flat rings. It does not have effect in this version.

* Openbabel