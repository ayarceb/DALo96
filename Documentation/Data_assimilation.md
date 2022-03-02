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

<img src="https://render.githubusercontent.com/render/math?math=y \sim \mathcal{N}\left(H\cdot x^*,R\right)

where  <img src="https://render.githubusercontent.com/render/math?math=B \in \Re^{n\times n}"> is the background error covariance matrix,  <img src="https://render.githubusercontent.com/render/math?math=H \in \Re^{m\times n}"> is a linear operator that propagates the state space into the observation space, and  <img src="https://render.githubusercontent.com/render/math?math=R \in \Re^{m\times m}"> is the observation error covariance matrix. 