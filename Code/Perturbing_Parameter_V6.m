% Synthetic data assimilation with the stochastic version of the Lorenz 96 (Perturbing the forcing parameter)
% Code to compare EnKF, EnKF schur, EnKF MC, EnKS data assimilation techniques for the Lorenz96 
% * the calculate of the gain in the next way  K2 = (Binv+ H'*(R\H)) for the analysis
clear all;close all;clc
% Folder where is located the taylor diagram utilities:

addpath '/Users/ayarcebotero/Dropbox/2020/ENKS_MC_paper/EnKS-MC/EnKS-MC_new'
%% ====================== Parameters for the model propagation, True generation ==========================
num=40;         % Number of ensembles
rng(30)
% load T_r_1.mat    % Load the schur product matrix radio=1
Tsim=400;          % Simulation time
dt=0.01;          % Step length
n=40;             % State number
m=n;             % Number of observations
ste_inter=30;     % State to focus the visualizations 
%F=8;             % Forcing factor model  In this version of the code the
%Forcing parameter is perturbed in this code version
radi=2;           % Radio for schur product
r=3;              % Predecesors radio for the Modified Cholesky
sigma=0.1;        % Noise factor for the error covariance matrix terms
R=sigma^2*eye(m);   % Observation Covariance matrix 
H = eye(n,n);       % Observation operator matrix
H = H(randperm(n,m),:);
M=2; %Ventana de asimilacion
frequency=1; % Frequency of observations  
pfactorEnkf = 1;  %INFLATION FACTOR
pfactorEnkfS = 1.2;  %INFLATION FACTOR
pfactorEnkfMC = 1.2;  %INFLATION FACTOR
%% ====================== Generate real state=========================================
% Part to generate the ensemble perturbing the forcing parameter. The part for visualize the forcing
% factor change incremental
x0=1*randn(n,1);incr=2;FF=7;
[Xreal]=Lorenz_96(Tsim+M*frequency,dt,x0,FF);            % Create vector
%% ====================== Data assimilation part =====================================
cont=1;
% for num=4:1:40    commented part for experiment with the ensembles
cont=cont+1;
Xb=zeros(n,num,Tsim);       % Background ensemble initialization
X_chol=zeros(n,num,Tsim);       % Background ensemble initialization
X_S=zeros(n,num,Tsim);
meanxb=zeros(n,Tsim);       % Background ensemble mean initialization
Xb(:,:,1)=2*randn(n,num);      % Initial condition creation

Xa=Xb;X_chol=Xb;X_S_chol=Xb;  
X_schur=Xb;                  % Initialize the Schur scenario background matrix
X_S=Xb;

meanxa_EnKF=zeros(n,Tsim);           % initialization of the vector to put the analisis ENKF mean
meanxa_EnKF_CHOLESKY=zeros(n,Tsim);  % initialization of the vector to put the analisis ENKF_CHOLESKY mean
meanxa_EnKF_schur=zeros(n,Tsim);     % initialization of the vector to put the analisis ENKF mean
meanxa_EnKS=zeros(n,Tsim);           % initialization of the vector to put the analisis EnKS mean
meanxa_EnKS_chol=zeros(n,Tsim);      % initialization of the vector to put the analisis EnKS mean
B=zeros(n,n,Tsim);            % initialization of the covariance background error matrix for EnKF
BC=zeros(n,n,Tsim);           % initialization of the covariance background error matrix for Cholesky
BS=zeros(n,n,Tsim);            % initialization of the covariance background error matrix for Smoother
%% ============================== Part for the observation frequency =================
muestreo=frequency:frequency:Tsim;     % vector to put the observation to assimilate
 
Y=H*Xreal;     %  Sampling of the observation operator
% Y=H*y;     %  Sampling of the observation operator
select=1; % this is used to select between the enkf (1), enkfMC (2), enks (3) and enksMC (4)
fact=0.6;  % factor for the forcing factor perturbation
%% ============================== Scenario EnKS ======================================
% if select==1
%   
 for i=2:Tsim-M*frequency
 %======== Forecast Step==========
 
    for en=1:num
      
      Dist=fact*randn(num,1);
      FF=FF+Dist(randi([1 num]),1); 
        [X_S(:,en,i)]=Lorenz_96_one_step(1,dt,squeeze(Xa(:,en,i-1)),FF);
    end
  
     meanxb(:,i)=mean(X_S(:,:,i),2);    
    
  
     XXb(:,:,i)= X_S(:,:,i)-meanxb(:,i);  
 
  


   B(:,:,i)=((1/(num-1))*(XXb(:,:,i)*XXb(:,:,i)'));


     
 


%     imagesc(XXb(:,:,i)*XXb(:,:,i)');colorbar;xlabel('states');ylabel('states')
% pause(1)
   if sum(muestreo==i)   %Validation if there are a observation available
%======================== Analysis Step ENKS=========================================

    X_aw(:,:,1)=X_S(:,:,i);
    
    if M==1
    B_aw(:,:,1)=B(:,:,i);    
    end
        for kk=2:(M*frequency)
            for en=1:num    
                [X_aw(:,en,kk)]=Lorenz_96_one_step(1,dt,squeeze(X_aw(:,en,kk-1)),FF);
            end
            meanxb_aw(:,kk)=mean(X_aw(:,:,kk),2);    
    
             for Nen=1:num
                XXb_aw(:,Nen,kk)= X_aw(:,Nen,kk)-meanxb_aw(:,kk)+1*randn(n,1);  
             end
  
             B_aw(:,:,kk)=((1/(num-1))*(XXb_aw(:,:,kk)*XXb_aw(:,:,kk)'));   
        end
    aux=0;
    for kk2=1:frequency:M*frequency
        K(:,:,kk2)=B_aw(:,:,kk2)*H'*pinv(H*B_aw(:,:,kk2)*H'+R);     % Matriz de Ganancia de Kalman
        C(:,kk2)= Y(:,i+kk2-1)+sigma*randn(m,1);
        aux=aux+(K(:,:,kk2)*(C(:,kk2)*ones(1,num)-H*X_aw(:,:,kk2)))/M;
    end
    
    for en=1:num
      Xa(:,en,i)=X_S(:,en,i)+aux(:,en);
    end
    

else 
  Xa(:,:,i)=X_S(:,:,i);  
end
  meanxa_EnKS(:,i)=mean(Xa(:,:,i),2);                                                                                                                                                                                                                                                                                                                                              
end
%% ============================== Scenario EnKS-MC ===================================
% % % %% 
% % % % if select==1
% % % %   
%  for i=2:Tsim-M*frequency
% % %======== Forecast Step==========
%  
% %      for en=1:num
% %           [X_S_chol(:,en,i)]=Lorenz_96_one_step(1,dt,squeeze(X_S_chol(:,en,i-1)),F);
% %      end
%  for en=1:num
% Dist=fact*randn(num,1);
%       FF=FF+Dist(randi([1 num]),1); 
%         [X_S_chol(:,en,i)]=Lorenz_96_one_step(1,dt,squeeze(Xa(:,en,i-1)),FF);
%  end
%         
%       meanxb(:,i)=mean(X_S_chol(:,:,i),2);    
%      
%    for Nen=1:num
%       XXb(:,Nen,i)= X_S_chol(:,Nen,i)-meanxb(:,i);  
%    end
%    
%  %  Standard Modified cholesky Background matrix 
% % 
%      % BC(:,:,i)=Calculo_B_Cholesky(XXb,r); %Estimation of Covariance by Modified Cholesky
%       BC(:,:,i)=(Calculo_B_Cholesky(XXb,r)); %Estimation of Covariance by Modified Cholesky
%       
% 
%  if sum(muestreo==i)   %Validation if there are a observation available
% % %======================== Analysis Step ENKS=========================================
% % 
%      X_aw(:,:,1)=X_S_chol(:,:,i);
%      
%      if M==1
%      B_aw(:,:,1)=BC(:,:,i);    
%      end
%          for kk=2:(M*frequency)
%              for en=1:num    
%                  Dist=fact*randn(num,1);
%       FF=FF+Dist(randi([1 num]),1); 
%         [X_aw(:,en,kk)]=Lorenz_96_one_step(1,dt,squeeze(Xa(:,en,i-1)),FF);
%                %  [X_aw(:,en,kk)]=Lorenz_96_one_step(1,dt,squeeze(X_aw(:,en,kk-1)),F);
%              end
%              meanxb_aw(:,kk)=mean(X_aw(:,:,kk),2);    
%      
%               for Nen=1:num
%                  XXb_aw(:,Nen,kk)= X_aw(:,Nen,kk)-meanxb_aw(:,kk)+0.2*randn(n,1);  
%               end
%    
%                B_aw(:,:,kk)=(Calculo_B_Cholesky(XXb_aw,r));   
%          end
%      aux=0;
%      for kk2=1:frequency:M*frequency
%          K(:,:,kk2)=B_aw(:,:,kk2)*H'*pinv(H*B_aw(:,:,kk2)*H'+R);     % Matriz de Ganancia de Kalman
%          C(:,kk2)= Y(:,i+kk2-1)+sigma*randn(m,1);
%          aux=aux+(K(:,:,kk2)*(C(:,kk2)*ones(1,40)-H*X_aw(:,:,kk2)))/M;
%      end
% %     
%      for en=1:num
%        X_S_chol(:,en,i)=X_S_chol(:,en,i)+aux(:,en);
%      end
%      
% % 
%  else 
%    X_S_chol(:,:,i)=X_S_chol(:,:,i);  
%  end
%    meanxa_EnKS_chol(:,i)=mean(X_S_chol(:,:,i),2);                                                                                                                                                                                                                                                                                                                                              
%  end
% % % 
% %% ============================== Scenario EnKF ======================================
% %if select==1
% for i=2:Tsim-1
% %======== Forecast Step==========
% % foplot(movmean(Error_ENKF,window));hold on
% % plot(movmean(Error_ENKF_schur,window));hold on
% % plot(movmean(Error_ENKF_MC,window));hold on
% % plot(movmean(Error_ENKS,window));hold on
% % plot(movmean(Error_ENKS_chol,window));hold on
%   for en=1:num
%       % FF=7;
%      Dist=fact*randn(num,1); FF=FF+Dist(randi([1 num]),1);
%      [Xb(:,en,i)]=Lorenz_96_one_step(1,dt,squeeze(Xb(:,en,i-1)),FF);
%    end
%   
%      meanxb(:,i)=mean(Xb(:,:,i),2);    
%      
% if sum(muestreo==i)   %Validation if there are a observation available
% %   for Nen=1:N
% %   XXb(:,Nen,i)= Xb(:,Nen,i)-meanxb(:,i);  
% %   end
%     Xbp(:,:,i)=meanxb(:,i)*ones(1,num)+ pfactorEnkf*(Xb(:,:,i)-meanxb(:,i)*ones(1,num));
%     meanXbp(:,i) = mean(Xbp(:,:,i),2);   
%     XXb(:,:,i) =Xbp(:,:,i)-meanXbp(:,i);
% % Standard EnKF Background matrix 
%     Binv=pinv((1/(num-1))*(XXb(:,:,i)*XXb(:,:,i)'));
% %     imagesc(XXb(:,:,i)*XXb(:,:,i)');colorbar;xlabel('states');ylabel('states')
% %     title(sprintf('Ensemble number %i',num))
% %     [L,D, Dsquare]=Calculo_B_CholeskyInc(XXb,r);
% %     Binv = L'*D*L;
% %======================== Analysis Step ENKF=========================================
% muestra=0;     
% %K=B*H'*pinv(H*B*H'+R);     % Matriz de Ganancia de Kalman
%   K2 = (Binv+ H'*(R\H));
%     for en=1:num
%          C(:,i)= Y(:,i)+sigma*randn(m,1); %  Almacenamiento de los datos sint??ticos
%           D= (C(:,i)-H*Xb(:,en,i));   
%          Xa(:,en,i)=Xb(:,en,i)+((K2\H')*(R\D));
%     end
%     Xb(:,:,i)=Xa(:,:,i);
% else
%      Xb(:,:,i)=Xb(:,:,i);
% end
% meanxa_EnKF(:,i)=mean(Xb(:,:,i),2);
% end
% %% ============================== Scenario EnKF_Schur_Localization ===================
% % 
% % if select==1
% %   
%  for i=2:Tsim-1
% % %======== Forecast Step==========
% %   for en=1:num
% %       [X_schur(:,en,i)]=Lorenz_96_one_step(1,dt,squeeze(X_schur(:,en,i-1)),F);
% %   end
%  
%    for en=1:num
%       % FF=7;
%       Dist=fact*randn(num,1);
%       FF=FF+Dist(randi([1 num]),1);
%      [X_schur(:,en,i)]=Lorenz_96_one_step(1,dt,squeeze(Xb(:,en,i-1)),FF);
%    end
%    
%   
%      meanxb(:,i)=mean(X_schur(:,:,i),2);    
%      
% if sum(muestreo==i)   %Validation if there are a observation available
% %      
% %    for Nen=1:N
% %       XXb(:,Nen,i)= X_schur(:,Nen,i)-meanxb(:,i);  
% %    end
% %   i
%   Xbp(:,:,i)=meanxb(:,i)*ones(1,num)+ pfactorEnkfS*(X_schur(:,:,i)-meanxb(:,i)*ones(1,num));
%   meanXbp(:,i) = mean(Xbp(:,:,i),2);   
%   XXb(:,:,i) =Xbp(:,:,i)-meanXbp(:,i)*ones(1,num);
% 
% % %  Standard EnKF Background matrix 
% % 
%     B=((1/(num-1))*(XXb(:,:,i)*XXb(:,:,i)'));
%     L=localization_matrix(n,radi);
%     B=pinv(L.*B); 
% % %======================== Analysis Step ENKF=========================================
% % muestra=0;     
% %K=B*H'*pinv(H*B*H'+R);     % Matriz de Ganancia de Kalman
%  
%   K2 = (B+ H'*(R\H));
%      for en=1:num
% %          
%           C(:,i)= Y(:,i)+sigma*randn(m,1);        %  Almacenamiento de los datos sint??ticos
%           D= (C(:,i)-H*Xb(:,en,i)); 
%           Xa(:,en,i)=Xb(:,en,i)+(K2\H')*(R\D);
%  
%      end
%      X_schur(:,:,i)=Xa(:,:,i);
%    
% else
%      X_schur(:,:,i)=X_schur(:,:,i);
%  end
% 
%    meanxa_EnKF_schur(:,i)=mean(X_schur(:,:,i),2); 
%  end
% % % 
% % 
% %% ============================== Scenario EnKF_MC ================================== 
 


%if select==2
  
for i=2:Tsim-1
%================= Forecast Step==========
   for en=1:num
     [X_chol(:,en,i)]=Lorenz_96_one_step(1,dt,squeeze(X_chol(:,en,i-1)),FF);
   end
   

     meanxb(:,i)=mean(X_chol(:,:,i),2);    

%      
%    for Nen=1:N
%       XXb(:,Nen,i)= X_chol(:,Nen,i)-meanxb(:,i);  
%    end
  Xbp(:,:,i)=meanxb(:,i)*ones(1,num)+ pfactorEnkfMC*(X_chol(:,:,i)-meanxb(:,i)*ones(1,num));
  meanXbp(:,i) = mean(Xbp(:,:,i),2);   
  XXb(:,:,i) =Xbp(:,:,i)-meanXbp(:,i)*ones(1,num);
%   i
%    %XXa=XXb(:,:,i)*sqrt(1/(N-1));   
%     BC(:,:,i)=Calculo_B_Cholesky(XXb,r); %Estimation of Covariance by Modified Cholesky
    inB=Calculo_inB_Cholesky(XXb,r); %Estimation of invserse Covariance by Modified Cholesky
% %  Bsquare=B^(1/2);
% %  L=Xb(:,:,i)-meanxb;
% %     
% %%  Standard EnKF MC Background matrix 
% 
%   
%      
% %=============== Analysis Step ENKF MC========================================     
if sum(muestreo==i)   %Validation if there are a observation available    
% KC=BC(:,:,i)*H'*pinv(H*BC(:,:,i)*H'+R);     % Matriz de Ganancia de Kalman para cholesky

 K2C = (inB+ H'*(R\H));     
%      
      for en=1:num
%          
          C(:,i)= Y(:,i)+sigma*randn(m,1);        %  Almacenamiento de los datos sint??ticos
%           Xa_C(:,en,i)=X_chol(:,en,i)+KC*(C(:,i)-H*X_chol(:,en,i));
          D= (C(:,i)-H*X_chol(:,en,i)); 
          Xa_C(:,en,i)= X_chol(:,en,i) + (K2C\H')*(R\D);
% 
      end
%      
      X_chol(:,:,i)=Xa_C(:,:,i);
    
 
%      

else
     X_chol(:,:,i)=X_chol(:,:,i);
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
 meanxa_EnKF_CHOLESKY(:,i)=mean(X_chol(:,:,i),2);

end
% %% ============================== Plots  ============================================
% %%% Taylor Diagrams
% % 
% % for j=1:n
% % [STATS] = taylor_statistics(meanxa_EnKF(j,1:i-20),Xreal(j,1:i-20));
% %      sdev_ENKF(j,:)=STATS.sdev;
% %      crmsd_ENKF(j,:)=STATS.crmsd;
% %      ccoef_ENKF(j,:)=STATS.ccoef;
% %      [STATS] = taylor_statistics(meanxa_EnKF_schur(j,1:i-20),Xreal(j,1:i-20));
% %      sdev_wrf(j,:)=STATS.sdev;
% %      crmsd_wrf(j,:)=STATS.crmsd;
% %      ccoef_wrf(j,:)=STATS.ccoef;
% %       [STATS] = taylor_statistics(meanxa_EnKF_CHOLESKY(j,1:i-20),Xreal(j,1:i-20));
% %      sdev_EnKF_MC(j,:)=STATS.sdev;
% %      crmsd_EnKF_MC(j,:)=STATS.crmsd;
% %      ccoef_EnKF_MC(j,:)=STATS.ccoef;
% %      
%       
% %      [STATS] = taylor_statistics(meanxa_EnKS(j,1:i-20),Xreal(j,1:i-20));
% %      sdev_EnKS(j,:)=STATS.sdev;
% %      crmsd_EnKS(j,:)=STATS.crmsd;
% %      ccoef_EnKS(j,:)=STATS.ccoef;
% %      [STATS] = taylor_statistics(meanxa_EnKS_chol(j,1:i-20),Xreal(j,1:i-20));
% %      sdev_EnKS_MC(j,:)=STATS.sdev;
% %      crmsd_EnKS_MC(j,:)=STATS.crmsd;
% %      ccoef_EnKS_MC(j,:)=STATS.ccoef;
% % end
% % label = containers.Map({'ENKF', 'ENKF Schur','ENKF MC','ENKS','ENKS MC'}, {'r.', 'b.','y.','c.','m.'});
% 
% % taylor_diagram(sdev_ENKF(:,2),crmsd_ENKF(:,2),ccoef_ENKF(:,2),...
% %     'markerLabel',label,  'markerKey', 'ENKF', ...
% %     'markerSize',15,...
% %     'tickRMS',0.0:2.0:10.0,'tickRMSangle',110.0, ...
% %     'colRMS','g', 'styleRMS', ':', 'widthRMS', 2.0, ...
% %     'tickSTD',0.0:1.0:4.0, 'limSTD',8.0, ...
% %     'colSTD',[127/255 0 1], 'styleSTD', '-.', 'widthSTD', 1.0, ...
% %     'colCOR','k', 'styleCOR', '--', 'widthCOR', 1.0);
% % figure
% % 
% % taylor_diagram(sdev_wrf(:,2),crmsd_wrf(:,2),ccoef_wrf(:,2),...
% %     'markerLabel',label,  'markerKey', 'ENKF Schur', ...
% %     'markerSize',15,...
% %     'tickRMS',0.0:2.0:10.0,'tickRMSangle',110.0, ...
% %     'colRMS','g', 'styleRMS', ':', 'widthRMS', 2.0, ...
% %     'tickSTD',0.0:1.0:4.0, 'limSTD',8.0, ...
% %     'colSTD',[127/255 0 1], 'styleSTD', '-.', 'widthSTD', 1.0, ...
% %     'colCOR','k', 'styleCOR', '--', 'widthCOR', 1.0);
% % 
% % figure
% % 
% % taylor_diagram(sdev_ENKF(:,2),crmsd_ENKF(:,2),ccoef_ENKF(:,2),...
% %     'markerLabel',label,  'markerKey', 'ENKF', ...
% %     'markerSize',15,...
% %     'tickRMS',0.0:2.0:10.0,'tickRMSangle',50, ...
% %     'colRMS','g', 'styleRMS', ':', 'widthRMS', 2.0, ...
% %     'tickSTD',0.0:1.0:4.0, 'limSTD',8.0, ...
% %     'colSTD',[127/255 0 1], 'styleSTD', '-.', 'widthSTD', 1.0, ...
% %     'colCOR','k', 'styleCOR', '--', 'widthCOR', 1.0);
% % 
% % 
% % taylor_diagram(sdev_wrf(:,2),crmsd_wrf(:,2),ccoef_wrf(:,2),...
% %      'markerLabel',label,  'markerKey', 'ENKF Schur', ...
% %     'markerSize',15,...
% %     'tickRMS',0.0:2.0:10.0,'tickRMSangle',50, ...
% %     'colRMS','g', 'styleRMS', ':', 'widthRMS', 2.0, ...
% %     'tickSTD',0.0:1.0:4.0, 'limSTD',8.0, ...
% %     'colSTD',[127/255 0 1], 'styleSTD', '-.', 'widthSTD', 1.0, ...
% %     'colCOR','k', 'styleCOR', '--', 'widthCOR', 1.0);
% % 
% % 
% % 
% % % taylor_diagram(sdev_wrf(:,2),crmsd_wrf(:,2),ccoef_wrf(:,2),...
% % %     'overlay','on','markerSize',15,'markerLabelColor', 'b', ...
% % %     'markerLabel',label, 'markerKey', 'ENKF Schur');
% % % 
% % 
% 
% % taylor_diagram(sdev_EnKF_MC(:,2),crmsd_EnKF_MC(:,2),ccoef_EnKF_MC(:,2),...
% %     'overlay','on','markerSize',15,'markerLabelColor', 'b', ...
% %     'markerLabel',label, 'markerKey', 'ENKF MC');
% % 
% % taylor_diagram(sdev_EnKS_MC(:,2),crmsd_EnKS(:,2),ccoef_EnKS(:,2),...
% %     'overlay','on','markerSize',15,'markerLabelColor', 'b', ...
% %     'markerLabel',label, 'markerKey', 'ENKS');
% % taylor_diagram(sdev_EnKS_MC(:,2),crmsd_EnKS_MC(:,2),ccoef_EnKS_MC(:,2),...
% %     'overlay','on','markerSize',15,'markerLabelColor', 'b', ...
% %     'markerLabel',label, 'markerKey', 'ENKS MC');
% %----------------------------------------------------------------------------------------------------------------
%%  Time propagation plot
window=10;Fx=15;Ft=15; % Fontsize Xlabel Ylabel 
% figure
set(gcf,'defaultTextInterpreter','latex');
subplot(2,4,1)
plot(Xreal(10,:),'Color','r','LineWidth',2)
hold on
% cont;
% plot(meanxa_EnKF(10,:),'--b','LineWidth',2);
% hold on
% plot(meanxa_EnKF_schur(10,1:i-20),'c','LineWidth',2);
% hold on
 plot(meanxa_EnKF_CHOLESKY(10,1:i-20),'-.g','LineWidth',2)
 hold on
plot(meanxa_EnKS(10,1:i-20),':y','LineWidth',3)
% % hold on
%  plot(meanxa_EnKS_chol(10,1:i-20),'-.k','LineWidth',2)
%  hold on
%
legend({'X truth','Xa EnKF MChol','Xa ENKS'})
% legend({'X truth','Xa ENKF','Xa ENKF Schur','Xa ENKS'})
% legend({'X truth','Xa ENKF'});
% legend({'X truth','Xa ENKF','xA ENKS'});
% legend({'X truth','Xa ENKF','Xa ENKS','Xa ENKS MChol'})
% legend({'X truth','Xa ENKF','Xa ENKF Schur','Xa ENKS','Xa ENKS MChol'})
%legend({'X truth','Xa ENKF','Xa ENKF Schur','EnKF MC','Xa ENKS','Xa ENKS MChol'})
title('State 10','FontSize',Ft,'interpreter','latex','interpreter','latex')
xlim([1 Tsim-window])
subplot(2,4,5)
for i=1:Tsim
% Error_ENKF(i)=norm(abs(sum(meanxa_EnKF(10,i)-Xreal(10,i))));hold on
% Error_ENKF_schur(i)=norm(abs(sum(meanxa_EnKF_schur(10,i)-Xreal(10,i))));hold on 
Error_ENKS(i)=norm(abs(sum(meanxa_EnKS(10,i)-Xreal(10,i))));hold on
% Error_ENKS_chol(i)=norm(abs(sum(meanxa_EnKS_chol(10,i)-Xreal(10,i))));hold on
 Error_ENKF_MC(i)=norm(abs(sum(meanxa_EnKF_CHOLESKY(10,i)-Xreal(10,i))));hold on
% set(gca, 'YScale', 'log')
end
%  plot(movmean(Error_ENKF,window),'Color',[0.5,0.02*cont,0.2]);hold on
%   plot(movmean(Error_ENKF,window));hold on
%   plot(movmean(Error_ENKF_schur,window));hold on
 plot(movmean(Error_ENKF_MC,window));hold on
 plot(movmean(Error_ENKS,window));hold on
%  plot(movmean(Error_ENKS_chol,window));hold on

%legend({'EnKF'})
%legend({'EnKF','EnKS'})
%legend({'EnKF','ENKF Schur','EnKS'})
%legend({'EnKF','EnKF schur ','EnKS','EnKS MChol'})
%legend({'EnKF','EnKF schur','EnKF MC','EnKS','EnKS MChol'})
xlim([1 Tsim-window]);title('norm of the residual')
set(gca, 'YScale', 'log')
subplot(2,4,2)
plot(Xreal(20,:),'r','LineWidth',2)
hold on
plot(meanxa_EnKF(20,:),'--b','LineWidth',2)
hold on
% plot(meanxa_EnKF_CHOLESKY(20,1:i-20),'-.g','LineWidth',2)
% hold on
plot(meanxa_EnKF_schur(20,1:i-20),'c','LineWidth',2)
 hold on
 plot(meanxa_EnKS(20,1:i-20),':y','LineWidth',3)
hold on
 plot(meanxa_EnKS_chol(20,1:i-20),':k','LineWidth',2)
% legend({'X truth','Xa ENKF','Xa ENKF Schur','Xa ENKS'})
% legend({'X truth','Xa ENKF'});
% legend({'X truth','Xa ENKF','Xa ENKS'})
% legend({'X truth','Xa ENKF','Xa ENKF Schur','Xa ENKS'})
 legend({'X truth','Xa ENKF','Xa ENKF Schur','Xa ENKS','Xa ENKS MChol'})
%legend({'X truth','Xa ENKF','Xa ENKF Schur','EnKF MC','Xa ENKS','Xa ENKS MChol'})
title('State 20','FontSize',Ft,'interpreter','latex')

xlim([1 Tsim-window])
subplot(2,4,6)
for i=1:Tsim
Error_ENKF(i)=norm(abs(sum(meanxa_EnKF(20,i)-Xreal(20,i))));hold on
Error_ENKF_schur(i)=norm(abs(sum(meanxa_EnKF_schur(20,i)-Xreal(20,i))));hold on
Error_ENKS(i)=norm(abs(sum(meanxa_EnKS(20,i)-Xreal(20,i))));hold on 
Error_ENKS_chol(i)=norm(abs(sum(meanxa_EnKS_chol(20,i)-Xreal(20,i))));hold on
% Error_ENKF_MC(i)=norm(abs(sum(meanxa_EnKF_CHOLESKY(20,i)-Xreal(20,i))));hold on
set(gca, 'YScale', 'log')
end
plot(movmean(Error_ENKF,window));hold on
plot(movmean(Error_ENKF_schur,window));hold on
% plot(movmean(Error_ENKF_MC,window));hold on
plot(movmean(Error_ENKS,window));hold on
plot(movmean(Error_ENKS_chol,window));hold on

%legend({'EnKF'})
%legend({'EnKF','EnKS'})
%legend({'EnKF','EnKS Schur','EnKS'})
legend({'EnKF','EnKF schur ','EnKS','EnKS MChol'})
%legend({'EnKF','EnKF schur','EnKF MC','EnKS','EnKS MChol'})
xlim([1 Tsim-window]);title('norm of the residual')

subplot(2,4,3)
plot(Xreal(30,:),'r','LineWidth',2);hold on
plot(meanxa_EnKF(30,:),'--b','LineWidth',2)
hold on
% plot(meanxa_EnKF_CHOLESKY(30,1:i-20),'-.g','LineWidth',2)
% hold on
plot(meanxa_EnKF_schur(30,1:i-20),'c','LineWidth',2)
hold on
plot(meanxa_EnKS(30,1:i-20),':y','LineWidth',3)
hold on
plot(meanxa_EnKS_chol(30,1:i-20),':k','LineWidth',2)
%legend({'X truth','Xa ENKF','Xa ENKF Schur','Xa ENKS'})
% %legend({'X truth','Xa ENKF','Xa ENKS'})
legend({'X truth','Xa ENKF','Xa ENKF Schur','Xa ENKS','Xa ENKS MChol'})
%legend({'X truth','Xa ENKF','Xa ENKF Schur','EnKF MC','Xa ENKS','Xa ENKS MChol'})

title('State 30','FontSize',Ft,'interpreter','latex');xlim([1 Tsim-window])
subplot(2,4,7)
for i=1:Tsim
Error_ENKF(i)=norm(abs(sum(meanxa_EnKF(30,i)-Xreal(30,i))));hold on
Error_ENKS(i)=norm(abs(sum(meanxa_EnKS(30,i)-Xreal(30,i))));hold on
Error_ENKF_schur(i)=norm(abs(sum(meanxa_EnKF_schur(30,i)-Xreal(30,i))));hold on
Error_ENKS_chol(i)=norm(abs(sum(meanxa_EnKS_chol(30,i)-Xreal(30,i))));hold on
% Error_ENKF_MC(i)=norm(abs(sum(meanxa_EnKF_CHOLESKY(30,i)-Xreal(30,i))));hold on
set(gca, 'YScale', 'log')
end
plot(movmean(Error_ENKF,window));hold on
plot(movmean(Error_ENKF_schur,window));hold on
%plot(movmean(Error_ENKF_MC,window));hold on
plot(movmean(Error_ENKS,window));hold on
plot(movmean(Error_ENKS_chol,window));hold on

% legend({'EnKF'})
%legend({'EnKF','EnKS'})
% legend({'EnKF','EnKS Schur','EnKS'})
legend({'EnKF','EnKF schur','EnKS','EnKS MChol'})
%legend({'EnKF','EnKF schur','EnKF MC','EnKS','EnKS MChol'})
xlim([1 Tsim-window]);title('norm of the residual')

subplot(2,4,4)
plot(Xreal(40,:),'r','LineWidth',2);hold on
plot(meanxa_EnKF(40,:),'--b','LineWidth',2)
hold on
plot(meanxa_EnKF_schur(40,1:i-20),'c','LineWidth',2)
% hold on
% plot(meanxa_EnKF_CHOLESKY(40,1:i-20),'-.g','LineWidth',2)
% hold on
  plot(meanxa_EnKS(40,1:i-20),':y','LineWidth',3)
% hold on
 plot(meanxa_EnKS_chol(40,1:i-20),':k','LineWidth',2)
% legend({'X truth','Xa ENKF','Xa ENKF Schur','Xa ENKS'})
%legend({'X truth','Xa ENKF','Xa ENKS'})
%legend({'X truth','Xa ENKF'})
legend({'X truth','Xa ENKF','Xa ENKF Schur','Xa ENKS','Xa ENKS MChol'})
%legend({'X truth','Xa ENKF','Xa ENKF Schur','EnKF MC','Xa ENKS','Xa ENKS MChol'})

title('State 40','FontSize',Ft,'interpreter','latex');xlim([1 Tsim-window])
subplot(2,4,8)
for i=1:Tsim
Error_ENKF(i)=norm(abs(sum(meanxa_EnKF(40,i)-Xreal(40,i))));hold on    
Error_ENKS(i)=norm(abs(sum(meanxa_EnKS(40,i)-Xreal(40,i))));hold on
Error_ENKF_schur(i)=norm(abs(sum(meanxa_EnKF_schur(40,i)-Xreal(40,i))));hold on
Error_ENKS_chol(i)=norm(abs(sum(meanxa_EnKS_chol(40,i)-Xreal(40,i))));hold on
%Error_ENKF_MC(i)=norm(abs(sum(meanxa_EnKF_CHOLESKY(40,i)-Xreal(40,i))));hold on
set(gca, 'YScale', 'log')
end
 plot(movmean(Error_ENKF,M));hold on
plot(movmean(Error_ENKF_schur,window));hold on
% plot(movmean(Error_ENKF_MC,window));hold on
 plot(movmean(Error_ENKS,window));hold on
plot(movmean(Error_ENKS_chol,window));hold on

% legend({'EnKF','EnKS'})
%legend({'EnKF','EnKS'})
% legend({'EnKF','EnKF Schur','EnKS'})
legend({'EnKF','EnKF schur','EnKS','EnKS MChol'})
%legend({'EnKF','EnKF schur','EnKF MC','EnKS','EnKS MChol'})
xlim([1 Tsim-window]);title('norm of the residual')

clear Xbp meanxb XXb
% end
%% ====================== Hovvmoller diagrams ========================================
% subplot(2,4,1)
% legend();
% figure
% 
% figure(3)
% subplot(1,3,1)
% imagesc((meanxa_EnKF(:,:)))
% colorbar
% subplot(1,3,2)
% % % imagesc(Xreal(:,:))
% colorbar