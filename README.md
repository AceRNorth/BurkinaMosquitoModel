# BurkinaMosquitoModel
For simulating mosquito populations in Burkina Faso
About:
The scripts in this repository were written to study genetic (gene-drive) methods of controlling mosquito vectors of malaria. The files are as follows.
DsxModel.cpp is c++ code for a simulation model of mosquito metapopulation dynamics. It is accompanied by the header file headDsxModel.h. The files RainData.csv and SetAll1.csv are input files that can be used to simulate the metapopulation model in a 10^6 km^2 area of West Africa containing Burkina Faso, represeting (historical) rain data, and settlement data respectively. The files NonSpatial_r2.wl and NonSpatial_r1_r2alleles.wl are Mathematica (Wolfram Research) scripts to simulate simpler versions of the model (without spatial structure); the first (NonSpatial_r2.wl) includes the possibility of non-functional alleles forming that are resistant to the gene drive, while the second also includes functional reistant alleles.
The c++ simulation codes have been tested with the Intel compiler v14.0.2 ; The Mathematica files have been tested in Mathematica v12.0.1
