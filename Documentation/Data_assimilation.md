# Data Assimilation (DA) 

is a mathematical process that is used to incorporate observations in a dynamic model,
to improve its representation of the reality.


The benefit of sequential DA methods is to estimate the current state <img src="https://render.githubusercontent.com/render/math?math=\mathbf{x}^* \in\Re^{n\times 1}"> of a dynamical system that evolves according to some numerical model operator, where <img src="https://render.githubusercontent.com/render/math?math=n"> is the number of states \cite{Evensen1994,Anderson1999},


   <img src="https://render.githubusercontent.com/render/math?math=\mathbf{x}^*_{k}=\mathcal{M}_{(k-1) \rightarrow k}( \mathbf{x}^{*}_{(k-1)})">


Where  <img src="https://render.githubusercontent.com/render/math?math=k"> denotes de time index, and the  <img src="https://render.githubusercontent.com/render/math?math=\mathcal{M}"> represents the model operator of the dynamics for instance. This estimation is performed based of a first guess or prior estimate  <img src="https://render.githubusercontent.com/render/math?math=x^b \in \Re^{n\times 1}$ of $\mathbf{x}^*">,
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

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{x}^a=\mathbf{x}^b%2B[\mathbf{P}^{-1}%2B\mathbf{H}^{T}\mathbf{R}^{-1}\mathbf{H}]^{-1}\mathbf{H}^{T}\mathbf{R}\mathbf{d}">


## EnKF Schur Product covariance localization

The emergence of misleading or spurious correlations between elements of the state space is inescapable due to the approximation of the state space covariance by a finite number of ensemble members. These spurious correlations can be removed by a precedura called localization. The covariance localization, also known as Schur localization, concentrates on the forecast error covariance matrix, removing longer-range correlations in the error covariances at a given distance. The pointwisemultiplication is called a Schur product and denoted by <img src="https://render.githubusercontent.com/render/math?math=\circ">:

<img src="https://render.githubusercontent.com/render/math?math=f   \circ    \boldsymbol{P}^f]_{i,j}=[\boldsymbol{P}^f]_{i,j}[f]_{i,j}">


## EnKF Modified Cholesky
## EnKS

The expression for the analysis given <img src="https://render.githubusercontent.com/render/math?math=n">  observations in the Data assimilation window for the smoother technique is given by the following expression


<img src="https://render.githubusercontent.com/render/math?math=\mathbf{x}^a=\mathbf{x}^b%2B\sum_{k=1}^s[\mathbf{B}_{0,k}^{-1}%2B\mathbf{H}_{k}^{T}\mathbf{R_k}^{-1}\mathbf{H}_{k}]^{-1}\mathbf{H}^{T}_{k}\mathbf{R_k}\mathbf{d}_{k}">


## EnKS Modified Cholesky


