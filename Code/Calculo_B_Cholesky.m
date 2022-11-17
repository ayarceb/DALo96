function [B,Bsquare]=Calculo_B_Cholesky(XB,r)
n=size(XB,1); %Number of states
T = eye(n,n);
D(1,1) = var(XB(1,:));
for i=2:n
   Neigh = mod([i-r:i-1 i+1:i+r],n); % Neighborhood according with the radius
   P = Neigh(Neigh<i); %Predecessors within the local domain
   P = P(P~=0);
   Xi=XB(i,:)';  %Current Ensemble member
   Xj=XB(P,:)';  % Predecessors ensemble members
%    Beta=-(pinv(Xj)*(Xi));  %Solution of Xi=sum_{j E P(i)} Xj*Betaj
   Beta=-(pinv(Xj'*Xj)*(Xj'*Xi));
   T(i,P) = Beta;  
   D(i,i) = var(Xi-Xj*Beta);
    
end
B = pinv(T)*D*pinv(T)';
Dsquare=sqrtm(D);
Bsquare=pinv(T)*Dsquare;

end