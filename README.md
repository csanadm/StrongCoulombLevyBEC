# StrongCoulombLevyBEC

## Basic info
This code is meant to calculate 2-pion Bose-Einstein correlation functions including the Coulomb and strong final state effects. This repository contains these directories and files:
- classes: base classes for the calculations
  - HypSInt_reader: 
  - Levy_reader: 
  - CorrFunc_reader:
  - CorrFunc_evaluator:
  - StrongCorrFunc_evaluator:
  - StrongCorrFunc_reader:
- sources: main calculation codes
  - HypSInt_filler: 
  - Levy_filler: 
  - CorrFunc_filler: 
  - functions: 
  - AfterBurner_CorrFunc_bettereta: 
  - plotter: 
- makefile: makefile for building the code
- Glue.sh: script to glue partial results together
- condor_inter.sh: running the calculations on condor

## Building
one just has to issue "make"

In addition, these directories are needed (to be calculated before doing "make"):
- tables: to store resulting data tables
- Data: to store resulting data
- Logs: to store logs
- Temp: to store temporary files
- object: to store object files
- object/deps: to store dependencies
- exe: to store final executables


## Running
TBA

