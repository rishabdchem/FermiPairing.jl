# FermiPairing.jl

Traditional and novel quantum chemistry methods for approximating exact energies of the 
reduced Bardeen--Cooper--Schrieffer (BCS) Hamiltonian, a model Hamiltonian describing fermion pairs.
Can be extended to any Hamiltonian made from operators that follow 
[su(2)](https://en.wikipedia.org/wiki/Lie_algebra) commutation relations.

## Methods

- Doubly occupied configuration interaction (DOCI).
- Particle-hole configuration interaction (CI) based on Hartree--Fock (HF).
- Antisymmetrized geminal power (AGP).
- Linear combinations of AGP (LC-AGP) as non-orthogonal configuration interaction (NOCI). 
- LC-AGP as non-orthogonal multi-configuration self-consistent field (NO-MCSCF). 
- Binary tree state (BTS). 
- Linear combinations of BTS (LC-BTS) as NOCI. 
- LC-BTS as NO-MCSCF. 

## Install 

There are two ways of adding Julia [packages](https://pkgdocs.julialang.org/v1/managing-packages/), 
either using the `add` command or the `dev` command. 
Use the Julia REPL, enter Pkg by pressing `]` and write the command below. 

```
pkg> add "git@github.com:rishabdchem/FermiPairing.jl.git"
```

Or clone this repo to your local machine, go to the Julia Pkg REPL, and write the command below.

```
pkg> dev path/to/FermiPairing.jl
```

## References

1. Dukelsky, Pittel, Sierra.
   Colloquium: Exactly solvable Richardson-Gaudin models for many-body quantum systems. 
   *Rev. Mod. Phys.* **76**, 643 (2004)
   [article](https://doi.org/10.1103/RevModPhys.76.643)
   [preprint](https://arxiv.org/abs/nucl-th/0405011)
 
1. Khamoshi, Henderson, Scuseria. 
   Efficient evaluation of AGP reduced density matrices.
   *J. Chem. Phys.* **151**, 184103 (2019)
   [article](https://doi.org/10.1063/1.5127850)
   [preprint](https://arxiv.org/abs/1909.06345)
 
1. Dutta, Henderson, Scuseria. 
   Geminal replacement models based on AGP. 
   *J. Chem. Theory Comput.* **16**, 6358 (2020) 
   [article](https://doi.org/10.1021/acs.jctc.0c00807)
   [preprint](https://arxiv.org/abs/2008.00552)
 
1. Dutta, Chen, Henderson, Scuseria. 
   Construction of linearly independent non-orthogonal AGP states. 
   *J. Chem. Phys.* **154** 114112 (2021)  
   [article](https://doi.org/10.1063/5.0045006)
   [preprint](https://arxiv.org/abs/2101.08911) 

1. Dutta, Gao, Khamoshi, Henderson, Scuseria. 
   Correlated pair ansatz with a binary tree structure. 
   [preprint](https://arxiv.org/abs/2310.20076) 
 
