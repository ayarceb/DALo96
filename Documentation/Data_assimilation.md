# Data Assimilation (DA) 

Is a mathematical process that is used to incorporate observations in a dynamic model,
to improve its representation of the reality.


The benefit of sequential DA methods is to estimate the current state <img src="https://render.githubusercontent.com/render/math?math=\mathbf{x}^* \in\Re^{n\times 1}"> of a dynamical system that evolves according to some numerical model operator, where <img src="https://render.githubusercontent.com/render/math?math=n"> is the number of states \cite{Evensen1994,Anderson1999},


   <img src="https://render.githubusercontent.com/render/math?math=\mathbf{x}^*_{k}=\mathcal{M}_{(k-1) \rightarrow k}( \mathbf{x}^{*}_{(k-1)})">


Where  <img src="https://render.githubusercontent.com/render/math?math=k"> denotes de time index, and the  <img src="https://render.githubusercontent.com/render/math?math=\mathcal{M}"> represents the model operator of the dynamics for instance. This estimation is performed based of a first guess or prior estimate  <img src="https://render.githubusercontent.com/render/math?math=x^b \in \Re^{n\times 1}$ of $\mathbf{x}^*">,Cancel changes
In ensembled-based methods, an ensemble of <img src="https://render.githubusercontent.com/render/math?math=N"> model realizations.


  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{X}_k^b=[\mathbf{x}_k^{b[1]},\mathbf{x}_k^{b[2]},..,\mathbf{x}_k^{b[N]}] \in \mathbb{R}^{n\times N}">
  

wich are assumed to be normally distributed


  <img src="https://render.githubusercontent.com/render/math?math=x \sim \mathcal{N}(x^b,B)">

Where the mean is assumed to be <img src="https://render.githubusercontent.com/render/math?math=x^b">
and a covariance matrix <img src="https://render.githubusercontent.com/render/math?math=B \in \Re^{n\times n}">. The observations are also assumed normal distributed 

<img src="https://render.githubusercontent.com/render/math?math=y \sim \mathcal{N}\left(H\cdot x^*,R\right)">

where  <img src="https://render.githubusercontent.com/render/math?math=B \in \Re^{n\times n}"> is the background error covariance matrix,  <img src="https://render.githubusercontent.com/render/math?math=H \in \Re^{m\times n}"> is a linear operator that propagates the state space into the observation space, and  <img src="https://render.githubusercontent.com/render/math?math=R \in \Re^{m\times m}"> is the observation error covariance matrix. 


## EnKF

The expression for the EnKF analysis is:

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{x}^a">

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{x}^a=\mathbf{x}^b%2B\sum_{k=1}^s[\mathbf{B}_{0,k}^{-1}%2B\mathbf{H}_{k}^{T}\mathbf{R_k}^{-1}\mathbf{H}_{k}]^{-1}\mathbf{H}^{T}_{k}\mathbf{R_k}\mathbf{d}_{k}">

where  <img src="https://render.githubusercontent.com/render/math?math=\mathbf{d}"> is the innovation term (difference between the observations and the model in the observation sites)

## EnKF Schur Product covariance localization

The emergence of misleading or spurious correlations between elements of the state space is inescapable due to the approximation of the state space covariance by a finite number of ensemble members. These spurious correlations can be removed by a precedura called localization. The covariance localization, also known as Schur localization, concentrates on the forecast error covariance matrix, removing longer-range correlations in the error covariances at a given distance. The pointwisemultiplication is called a Schur product and denoted by <img src="https://render.githubusercontent.com/render/math?math=\circ">:


![App Lorenz 96](https://github.com/ayarceb/Data-Assimilation-interactive-tool/blob/main/Localization.png)

## EnKF Modified Cholesky


<img src="https://render.githubusercontent.com/render/math?math=\mathbf{x}^a=\mathbf{x}^b%2B\sum_{k=1}^s[\mathbf{B}_{0,k}^{-1}%2B\mathbf{H}_{k}^{T}\mathbf{R_k}^{-1}\mathbf{H}_{k}]^{-1}\mathbf{H}^{T}_{k}\mathbf{R_k}\mathbf{d}_{k}">


The computation of the inverse covariance matrix via MC has two fold benefits: the former is the contruction of a localized matrix from the coefficient calculations and the generation of a sparse matrix that for large scale problems may benefit calculation and storage.

![App Lorenz 96](https://github.com/ayarceb/Data-Assimilation-interactive-tool/blob/main/Cholesky_Dalo96.png)

![App Lorenz 96](https://github.com/ayarceb/Data-Assimilation-interactive-tool/blob/main/Cholesky_Dalo96_2.png)







## EnKS

The expression for the analysis given <img src="https://render.githubusercontent.com/render/math?math=n">  observations in the Data assimilation window for the smoother technique is given by the following expression



## EnKS Modified Cholesky


