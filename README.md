![kép](https://github.com/user-attachments/assets/3b0997d4-08c5-496e-877a-1fe4b72b25a1)# StrongCoulombLevyBEC

## Basic info
This code is meant to calculate 2-pion Bose-Einstein correlation functions including the Coulomb and strong final state effects, more specifically equation (37) from [arXiv:1912.01381](https://arxiv.org/abs/1912.01381). The calculation proceeds through creating binary tables and using them in subsequent steps of the calculation (in order to speed up the process). The code also utilizes condor as a job distribution framework, is hence computing environment specific.

This repository contains these directories and files:
- `sources`: main calculation codes
  - `functions`: definitions of special functions needed for the calculations
  - `Levy_filler`: calculates the D(r) source function according to equation (11) from [arXiv:1912.01381](https://arxiv.org/abs/1912.01381) in a distributed manner, creates `Levy_table_<N>.dat` files (where `<N>` represents the job number), which then have to be concatenated to `Levy_table.dat`
  - `HypSInt_filler`: calculates the dy integrals, creates `HypSInt_table_<N>.dat` files (again `<N>` represents the job number), which then have to be concatenated to `HypSInt_table.dat`
  - `CorrFunc_filler`: calculates the dr integrals using `HypSInt_table.dat` and `Levy_table.dat`, creates `CorrFunc_raw_table_<N>.dat` files, which then have to be concatenated to `CorrFunc_raw_table.dat`
  - `CorrFunc_filler_kaons`: same for kaons
  - `CorrFunc_filler_kaons`: same for protons
  - `AfterBurner_CorrFunc_bettereta`: Combines the results of the four integrals from equations (39)-(42) from arXiv:1912.01381 according to equation (38) using `CorrFunc_raw_table.dat` and creates `CorrFunc_table.dat`
  - `AfterBurner_CorrFunc_kaons`: same for kaons
  - `AfterBurner_CorrFunc_protons`: same for protons
  - `plotter`:  just a simple toolset to make plotting with root easier
- `classes`: base classes for reading the results of the calculations
  - `Levy_reader`: used to read the output of Levy_filler
  - `HypSInt_reader`: used to read the output of HypSInt_filler
  - `StrongCorrFunc_reader`: used to read the output of CorrFunc_filler
  - `StrongCorrFunc_evaluator`: used to evaluate the correlation function for any values of the source parameters and the realtive momentum variable by using interpolation.
- `makefile`: makefile for building the code
- `condor_inter.sh`: running the calculations on condor
- `Glue.sh`: script to glue (concatenate) partial results (of the individual condor jobs) together

## Building
One just has to issue "make"

In addition, these directories are needed (to be calculated before doing "make"):
- `tables`: to store resulting data tables
- `Data`: to store resulting data
- `Logs`: to store logs
- `Temp`: to store temporary files
- `object`: to store object files
- `object/deps`: to store dependencies
- `exe`: to store final executables


## Running
One has to first create the Levy tables via `Levy_filler` and the HypSInt tables via `HypSInt_filler`, then concatenate them via `Glue.sh` (creating a list in the specific order, for example using `ls -1v`). One then runs `CorrFunc_filler` to create CorrFunc tables, and concatenates those as well. Finally, one creates the final tables by running `AfterBurner_CorrFunc_bettereta` as well.


## Example plot
Below we show an example plot with correlation functions of identical charged pion, kaon and proton pairs, for R=3 fm and alpha=1
![exampfunc_PrPiKa](https://github.com/user-attachments/assets/d9b53aff-fddb-42ef-9da6-7271d887b6d9)
