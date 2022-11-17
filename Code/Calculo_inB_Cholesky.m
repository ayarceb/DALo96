function inB=Calculo_inB_Cholesky(XB,r)
n=size(XB,1); %Number of states
T = eye(n,n);
inD(1,1) = 1/var(XB(1,:));
for i=1:n
   Neigh = mod([i-r:i-1 i+1:i+r],n); % Neighborhood according with the radius
   P = Neigh(Neigh<i);  %Predecessors within the local domain
   P = P(P~=0);
   Xi=XB(i,:)';  %Current Ensemble member
   Xj=XB(P,:)';  % Predecessors ensemble members
   Beta=-(pinv(Xj'*Xj)*(Xj'*Xi));
   T(i,P) = Beta;
   inD(i,i) = 1/var(Xi-Xj*Beta);
   
end
inB = T'*inD*T;

end