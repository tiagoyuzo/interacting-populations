%--------------------------------------------------------------------------
%   parametros
%
%   define os parametros do problema:
%
%   alfa(i): coeficiente de difusao, i = 1,2,3
%   vx(i): velocidade em x, i = 1,2,3
%   vy(i): velocidade em y, i = 1,2,3
%   sigma: decaimento do poluente
%   lambda(i): taxa de reproducao intrinseca, i = 2,3
%   kp(i): capacidade de suporte, i = 2,3
%   rho(i): coeficiente de escala para relacao com c, i = 2,3
%   nuc(i): coeficiente de mortalidade causado por c, i = 2,3
%   nup(i): coeficiente de mortalidade causado por p ou q, i = 2,3
%   beta(i): coeficiente da condicao de contorno, i = 1,2,3
%
%   i = 1: parametros de c
%   i = 2: parametros de p
%   i = 3: parametros de q
%--------------------------------------------------------------------------
function [alfa,vx,vy,sigma,lambda,kp,rho,nuc,nup,beta] = parametros
%---======================================---------------------------------
%   parametros da concentracao de poluicao
%---======================================---------------------------------
alfa(1) = 0.2;
vx(1) = 0.2;
vy(1) = 0.2;
sigma = 4; %sigma = 4;
beta(1) = 0.2;
%---=========================----------------------------------------
%   parametros da populacao p
%---=========================----------------------------------------
alfa(2) = 0.002;
vx(2) = 0.002;
vy(2) = 0.002;
lambda(2) = 2;
kp(2) = 100;
rho(2) = 200;
nuc(2) = 100; %50;
nup(2) = 0.002;
beta(2) = 0.02;
%---=========================----------------------------------------
%   parametros da populacao q
%---=========================----------------------------------------
alfa(3) = 0.002;
vx(3) = 0.002;
vy(3) = 0.002;
lambda(3) = 5;
kp(3) = 100;
rho(3) = 200;
nuc(3) = 100; %50;
nup(3) = 0.002;
beta(3) = 0.02;
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------