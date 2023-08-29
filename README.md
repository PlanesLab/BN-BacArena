# BN-BacArena
Bayesian Network implementation in BacArena. This algorithm allows the dynamic modelling and simulation of bacterial communities considering the microbe-microbe and microbe-nutrient interactions. For more information, please refer to [paper]

The contents of BN-BacArena can be downloaded by cloning this repository.
```
git clone https://github.com/PlanesLab/BN-BacArena.git
```

## Contents
- **R** and **MATLAB** folders contain the source code required to launch BN-BacArena algorithm in both programming languages respectively.
- **test** folder contains the necessary data to launch BN-BacArena in both programming languages.
- The scripts **testR.R** and **testMat.m** allow to automatically run BN-BacArena algorithm in R and MATLAB respectively using data located at **test** folder.

## Requirements
Prior installation of software/library is required to run BN-BacArena and/or test scripts.
- R: Please install **BacArena** and **openxlsx** libraries prior running BN-BacArena.
- MATLAB: The installation of **The COBRA Toolbox** is required to run BN-BacArena. More information available [here](https://opencobra.github.io/cobratoolbox/latest/installation.html).

## Notes
- Note that random seeds are different between **R** and **MATLAB** languages. Therefore, BN-BacArena results could be different between them even if the same seed is applied.
- Using coefficient matrices is not a requirement for employing the BN-BacArena code (in both R and MATLAB). By doing so, users can still achieve the identical results provided by the original BacArena code. Accordingly, cell type and nutrient coefficient matrices can be applied independently.
- To use your own coefficient matrices, please remember that the dimensions of the matrices are ***c***x***c*** and ***c***x***n***, where ***c*** and ***n*** are the different cell types introduced in the arena and the nutrients of the culture media respectively.
- In order to use a nutrient coefficient matrix in **R**, it is mandatory to add the argument ```addAnyWay = TRUE``` in the **addSubs** function of BacArena.
```
arenaWithMedia = addSubs(arena, ..., addAnyWay = TRUE)
```
