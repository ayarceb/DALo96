# Data-Assimilation-interactive-tool for the Lorenz 96 
The interactive tool DALo96(V0.1) lets explore the Lorenz 96 dynamics changing the simulation configuration interactively.

The Lorenz 96 model is defined by

![img](http://latex.codecogs.com/svg.latex?%5Cfrac%7Bdx_i%7D%7Bdt%7D%3D%28x_%7Bi%2B1%7D-x_%7Bi-2%7D%29x_%7Bi-1%7D-x_i%2BF)

where the index  ![img](http://latex.codecogs.com/svg.latex?i) is cyclic which means ![img](http://latex.codecogs.com/svg.latex?x_%7B-1%7D%3Dx_%7Bn-1%7D), ![img](http://latex.codecogs.com/svg.latex?x_%7B0%7D%3Dx_%7Bn%7D)   and   ![img](http://latex.codecogs.com/svg.latex?x_%7Bn%2B1%7D%3Dx_%7B1%7D). 


DALo96 main purpose is to be a tool for explain general concepts of different Data Assimilation (DA) techniques using the Lorenz 96 model. ThCancel changese DA techniques that is possible to explore are

- EnKF
- EnKF Schur Product covariance localization
- EnKF Modified Cholesky
- EnKS
- EnKS Modified Cholesky

For the different implementation the forcing parameter <img src="https://render.githubusercontent.com/render/math?math=F"> is disturbed with an aditive noise
<img src="https://render.githubusercontent.com/render/math?math=\sim N(F_0,\gamma)"> to generate the ensemble space that promotes freedom degrees for a number of model propagations.
The forcing parameter <img src="https://render.githubusercontent.com/render/math?math=F_0"> can be adjusted in the slider and it can be disturbed with a multiplicative factor <img src="https://render.githubusercontent.com/render/math?math=\gamma"> with a knob.

![App Lorenz 96](https://github.com/ayarceb/Data-Assimilation-interactive-tool/blob/main/front.png)

App developed in App designer in Matlab R2021b
