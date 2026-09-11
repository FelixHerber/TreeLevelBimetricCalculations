# TreeLevelBimetricCalculations
Code for computing the tree-level scattering amplitudes of ghost-free bimetric theory.
Also contains an expansion of bimetric theory to quartic order in the mass eigenstates using xAct.
Uses Mathematica together with the packages xAct, FeynRules, FeynArts and FeynCalc.


## Getting started

Readers interested only in computing amplitudes can begin with
`ExampleCalculationsBimetricAmplitudes.nb`. 

Readers who only want to browse the amplitudes can instead use BimetricAmplitudeBrowser.nb.
The required model files and intermediate expansions have already been generated. 

Readers interested in reproducing the full pipeline can instead begin with
`xActExpansionOfBimetric.nb` and follow the workflow shown below.

## Description of the files

Workflow of the Notebooks:

``` xActExpansionOfBimetric.nb
          │
          ▼
xActExpansionOfBimetricTheory.fr
          │
          ▼
frForBimetric.nb
          │
          ▼
FeynArts model
          │
          ▼
BimetricAmplitudes.wl
          │
          ▼
Example notebook / HPC computation
```

### xActExpansionOfBimetric.nb
Contains the complete expansion of bimetric theory around flat backgrounds to fourth order in the perturbations. The expansion is given both in terms of the perturbations of g and f and in terms of the mass eigenstates G and M. All steps in obtaining this result, such as expanding the Einstein-Hilbert action terms or the bimetric interaction terms, are shown. At the end of the notebook, some code exists which converts the expansion of bimetric theory into text. The resulting text is written in a form that can be read by FeynRules and is the Lagrangian in the xActExpansionOfBimetricTheory.fr file.

Requires: 
- xAct and associated packages referenced below.

### xActExpansionOfBimetricTheory.fr
The FeynRules model file which contains the Lagrangian of bimetric theory expanded to fourth order in the mass eigenstates and written in a form FeynRules can understand.

### frForBimetric.nb
The notebook that uses FeynRules to read xActExpansionOfBimetricTheory.fr and then generates the FeynArts model files that can be found in the xActExpansionOfBimetricTheory_FA folder. Some checks are also applied in this notebook, such as checking whether the Lagrangian is Hermitian and whether it is diagonal in the mass eigenstates. The notebook contains absolute paths from the original development environment. Before running it on another system, replace the paths to the FeynRules installation and the .fr model file with the corresponding local paths.

This program generates the xActExpansionOfBimetricTheory_FA folder.

Requires: 
- FeynRules

Reads: 
- xActExpansionOfBimetricTheory.fr


### BimetricAmplitudes.wl
This contains all functions that allow for efficient computation of tree-level 2->2 helicity amplitudes in bimetric theory using FeynCalc.

This file was developed incrementally during the project. Some functions therefore use different scoping and return styles. The current version has been retained because it has been tested against the calculations used in the project.

Requires: 
- FeynCalc

### allBimetricAmplitudes.wl
Contains a list of all 199 distinct 2->2 helicity amplitudes of bimetric theory, where distinct means distinct with respect to the symmetries of time reversal, parity and initial- and final-particle exchange. See section 4.2 of https://arxiv.org/pdf/2609.10760 for an explanation

### ExampleCalculationsBimetricAmplitudes.nb
The example notebook which shows how to use the functions in BimetricAmplitudes.wl.
This code must be in the same folder as: BimetricAmplitudes.wl, allBimetricAmplitudes.wl and xActExpansionOfBimetricTheory_FA.
                                           
Requires: 
- FeynCalc
                                           
### BimetricAmplitudesHPC.wl 
Computes all 199 distinct tree-level helicity amplitudes using the functions in BimetricAmplitudes.wl on an HPC. 
This script is intended to be run on a cluster using the SLURM workload manager.

This code must be in the same folder as: BimetricAmplitudes.wl, xActExpansionOfBimetricTheory_FA and a folder SavedAmplitudes.
The SavedAmplitudes folder stores the individual amplitudes as .wl files.

Requires:
- FeynCalc

### GetOrbitsBinAmps.nb
This notebook generates representatives of the 199 distinct orbits of the group G. 

### BimetricAmplitudeBrowser.wl
This .wl contains functions allowing for browsing the amplitudes in allBimetricAmplitudes.wl.
For example, the functions allow amplitudes outside the set of 199 representatives to be reconstructed using the symmetry transformations.
Uses the transformations described in section 4.2 of our paper https://arxiv.org/pdf/2609.10760 to do this.

### BimetricAmplitudeBrowser.nb
This notebook is the user interface for BimetricAmplitudeBrowser.wl. 
The user can in this notebook browse through all possible 2->2 tree-level amplitudes of bimetric theory here, even those
that are not part of the 199 computed.

This code must be in the same folder as: BimetricAmplitudes.wl and BimetricAmplitudeBrowser.wl.

## Important note regarding the spin-2 propagator

The FeynArts folder generated by FeynRules will need to be modified to be able to use it.
The reason is the .gen file in this folder will contain this spin-2 propagator:
           
```Mathematica 

AnalyticalPropagator[Internal][ s1 T[j1, mom, {li1p1, li1p2} -> {li2p1, li2p2}] ] == (I*(((-((FourVector[mom, li1p1]*FourVector[mom, li2p2])/Mass[T[j1]]^2) + MetricTensor[li1p1, li2p2])*(-((FourVector[mom, li1p2]*FourVector[mom, li2p1])/Mass[T[j1]]^2) + MetricTensor[li1p2, li2p1]))/2 + ((-((FourVector[mom, li1p1]*FourVector[mom, li2p1])/Mass[T[j1]]^2)+ MetricTensor[li1p1, li2p1])*(-((FourVector[mom, li1p2]*FourVector[mom, li2p2])/Mass[T[j1]]^2) + MetricTensor[li1p2, li2p2]))/2 -((-((FourVector[mom, li1p1]*FourVector[mom, li1p2])/Mass[T[j1]]^2) + MetricTensor[li1p1, li1p2])*(-((FourVector[mom, li2p1]*FourVector[mom, li2p2])/Mass[T[j1]]^2) + MetricTensor[li2p1, li2p2]))/3)) PropagatorDenominator[mom, Mass[T[j1]]]

```
FeynRules generates the massive spin-2 propagator by default. This propagator cannot be used for the massless spin-2 field G, since the massless propagator is not obtained simply by taking the zero-mass limit of the massive spin-2 propagator, as reflected by the van Dam–Veltman–Zakharov discontinuity.

Instead, these lines must be replaced by
```Mathematica 
AnalyticalPropagator[Internal][ s1 T[j1, mom, {li1p1, li1p2} -> {li2p1, li2p2}] ] == FAPropagatorNumerator[li1p1, li1p2, li2p1, li2p2, mom, x]
```
Note that the model file: xActExpansionOfBimetricTheory_FA contains the correct propagator.

Special code must then be written to deal with FAPropagatorNumerator, which I have done in BimetricAmplitudes.wl.

## Tested versions

The repository has been tested with:

- Mathematica 14.3
- FeynCalc 10.1.0
- FeynArts 3.12
- FeynRules 2.3
- xPerm 1.2.3
- xTensor 1.2.0
- xPert 1.0.6
- TexAct 0.4.3
- xCoba 0.8.6

Later versions may also work, but have not been tested.
The files BimetricAmplitudeBrowser.nb and BimetricAmplitudeBrowser.wl were developed by Joakim Flinckman, while all other files were developed and maintained by Felix Herber, originating from his Master’s degree project.
Because Felix Herber's current research concerns other topics, active maintenance and support may be limited.

## References
This work would not have been possible without xAct, FeynRules, FeynArts, and FeynCalc.
Below are the references to the corresponding packages.

### xAct
This work uses the xAct tensor computer algebra packages developed by José M. Martín-García et al. 

Reference:
Nutma, T., & Mühlenberg, A. (2014). xTras: A field-theory inspired xAct package for Mathematica. Computer Physics Communications, 185, 1719-1738. https://doi.org/10.1016/j.cpc.2014.02.006 (arXiv:1308.3493)

Website: http://www.xact.es/

### FeynRules
Reference:
Alloul, A., et al. (2014). FeynRules 2.0 - A complete toolbox for tree-level phenomenology. Comput. Phys. Commun., 185, 2250-2300. DOI: 10.1016/j.cpc.2014.04.012

Website: https://cp3.irmp.ucl.ac.be/projects/feynrules/wiki for how to install it.

### FeynArts
Reference:
Hahn, T., "Generating Feynman diagrams and amplitudes with FeynArts 3," Comput. Phys. Commun. 140 (2001) 418-431, DOI: 10.1016/S0010-4655(01)00290-9.

Website: https://feynarts.de/ 

### FeynCalc
Reference:
Shtabovenko, V., Mertig, R., & Orellana, F. (2024). FeynCalc 10: Do multiloop integrals dream of computer codes? Computer Physics Communications, 109357. https://doi.org/10.1016/j.cpc.2024.109357

Website: https://feyncalc.github.io/

### The Master's Thesis

F. Herber (2025) The Scattering Amplitudes of Bimetric Theory. [Master's thesis, Stockholm University]. Diva Portal. 
Link: https://urn.kb.se/resolve?urn=urn:nbn:se:su:diva-248139

### arXiv
F. Herber, J. Flinckman, and S.F. Hassan (2026). Complete set of tree - level 2 -> 2 scattering amplitudes of ghost - free bimetric theory. arXiv. Link: https://arxiv.org/pdf/2609.10760

# Citation

This repository accompanies the paper below. If you use the code, results, or scattering amplitudes from this repository in academic work, please cite the corresponding paper:

F. Herber, J. Flinckman, and S. F. Hassan, Complete set of tree-level 2 -> 2 scattering amplitudes of ghost-free bimetric theory, arXiv:2609.10760 (2026).

https://arxiv.org/abs/2609.10760