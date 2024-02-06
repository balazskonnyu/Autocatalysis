[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10560157.svg)](https://doi.org/10.5281/zenodo.10560157)

# Autocatalysis

The following seven files contain all the program codes the results of which have been presented in [] . 

Files:
1. **_simualtion_controller.sh_** This file controls all types of simulations described in the paper. All the parameters are specified along with a detailed description.
2. **_parameter_generator_v3.c_**. This file generates the _parameters.txt_ file based on the instructions in the _simulaton_controller.sh_.

3. Parameters in the **_parameters.txt_** file 
* Section “general info”
    - start: start of simulation
    -  end: end of simulation
    - DeltaT: print the data at given time intervals (into res.res)
    - Num. of reactions: the number of reactions
    - Num. of Chem. Entities: number of chemical species
    - Sum of cc: total concentration of of chemical species
    - V: volume of jar in chemostate
    - v: in- and outflux rate in chemostate
    - Tmodel: type of simulation (0: resource-unlimited system or 1: chemostate)
    - JoinedCycle: cycles is joined (1) or not (0)
    - Num. of Interaction Type 1: number of monomolecular cross-catalysis
    - Num. of Interaction Type 2: number of bimolecular cross-catalysis
* Section “inflow”: list of the concentrations of the chemical species in the inflow, in the following form: \<type of chemical species\>\<comma\>\<inflow concentration\>
* Section “initial cc”: list of the initial concentrations of all chemical species of the system, in the following form: \<type of chemical species\>\<comma\>\<concentration\>
* Section “reactions”: list of the chemical reactions in the system, in the following form: \<Educt1\>\<comma\>\<Educt2\>\<comma\>\<k forward\>\<comma\>\<k backward\>\<comma\>\<Product1\>\<comma\>\<Product2\>
* Section “Interaction type 1”: list of monomolecular cross catalytic reaction types, in the following form: \<Substrate\>\<comma\>\<Enzyme\>\<comma\>\<kon\>\<comma\>\<koff\>\<comma\>\<kcat forward\>\<comma\>\<kcat backward\>\<comma\>\<Product1\>\<comma\>\<Product2q>
* Section “Interaction type 2”: list of bimolecular cross catalytic reaction types, in the following form: \<Substrate1\>\<comma\>Substrate2\>\<comma\>\<Enzyme\>\<comma\>\<kon1\>\<comma\>\<koff1\>\<comma\>\<kon2\>\<comma\>\<koff2\>\<comma\>\<kcat forward\>\<comma\>\<kcat backward\>\<comma\>\<Product1\>\<Product2\>
* Section “Joined Cycles at these points”: list of the common members of the cycles, in the following form: \<Member of CycleA=Member of CycleB\>

4. **_parsing_v3.c_**. This file parses the _parameters.txt_ file, creates the corresponding differential equation, and inserts it into the _compactAutCat.c_ file. The output of _compactAutCat.c_ has introduced in the scientific paper as the results of this project.
 
5. **_Kset_classification_None.ipynb_** contains the _PyTorch_ neural network codes in simulations of non-regulated cycles (_„None”_)

6. **_Kset_classification_Rev.ipynb_** contains the _PyTorch_ neural network codes in simulations with the last reaction step _reversible_.

7. **_Kset_classification_HetCat.ipynb_** contains the _PyTorch_ neural network codes in simulations with a member of one cycle catalysing a reaction of the other cycle (_cross-catalysis_).
 
