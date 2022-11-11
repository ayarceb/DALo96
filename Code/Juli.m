clc;close all;clear all

% Este codigo calcula un aproximado de los costos de traduccion

H=1;    % Horas de la pelicula
T=60*H;  % Tiempo total en minutos de la pelicula

% Tasa de tiempo real invertido 1 a 3 (por 1 minuto de audio de la pelicula cuanto dedico) 

Tasa=3;

TiempoReal=T*Tasa;   % Tiempo pelicula por tiempo invertido  

Thoras_real=TiempoReal/60;

DiasTrabajo=Thoras_real/8;   % DIAS de trabajo (HORAS LABORALES 8)

sprintf('Dias de trabajo %1.2f',DiasTrabajo)

% Costo dia de trabajo pesos (140.000 COP)

CT=220000;  % Poner aqui la variable de cuanto cobrarias por dia

CD=DiasTrabajo*CT;

sprintf('Ganancia proyecto %1.f',CD)