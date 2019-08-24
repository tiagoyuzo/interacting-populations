%--------------------------------------------------------------------------
%   submatriz - parte constante
%
%   constroi a parte linear das submatrizes e dos subvetores
%
%   programa baseado no de Pulino, P., obtido em:
%   http://www.ime.unicamp.br/~pulino/MT623/programas/
%--------------------------------------------------------------------------
function [meck,mdck,mepk,mdpk,meqk,mdqk,bck,bpk,bqk] = submatriz(M,C,k)
%--------------------------------------------------------------------------
global dt alfa vx vy sigma lambda
%--------------------------------------------------------------------------
%   declaracao das submatrizes e subvetores
%
%   meck: k-esima submatriz a esquerda de c
%   mdck: k-esima submatriz a direita de c
%   mepk: k-esima submatriz a esquerda de p
%   mdpk: k-esima submatriz a direita de p
%   meqk: k-esima submatriz a esquerda de q
%   mdqk: k-esima submatriz a direita de q
%   bck: k-esimo subvetor do lado direito de c
%   bpk: k-esimo subvetor do lado direito de p
%   bqk: k-esimo subvetor do lado direito de q
%--------------------------------------------------------------------------
meck = zeros(3,3);
mdck = zeros(3,3);
mepk = zeros(3,3);
mdpk = zeros(3,3);
meqk = zeros(3,3);
mdqk = zeros(3,3);
bck = zeros(3,1);
bpk = zeros(3,1);
bqk = zeros(3,1);
%--------------------------------------------------------------------------
%   pontos e pesos da quadratura para o triangulo padrao
%--------------------------------------------------------------------------
[w,t,nptos] = quadratura(1);
%--------------------------------------------------------------------------
%   parametros do modelo
%--------------------------------------------------------------------------
% [alfa,vx,vy,sigma,lambda,~,~,~,~,~] = parametros;
%--------------------------------------------------------------------------
%   calculo de constantes auxiliares
%--------------------------------------------------------------------------
dt2 = dt/2;
%--------------------------------------------------------------------------
%   loop de cada ponto de integracao
%--------------------------------------------------------------------------
for m = 1:nptos
    %----------------------------------------------------------------------
    %   integracao numerica: xi e eta sao os pontos de integracao
    %----------------------------------------------------------------------
    xi = t(1,m);
    eta = t(2,m);
    %----------------------------------------------------------------------
    %   funcoes de base padrao e seus gradientes
    %----------------------------------------------------------------------
    [phi,grad] = base(xi,eta);
    %----------------------------------------------------------------------
    %   transformacao afim: coordenadas reais x e y; (B^t)^-1 e 
    %                       jacobiano da transformacao
    %----------------------------------------------------------------------
    [Btinv,jacob] = transf(M,C,k);
    %----------------------------------------------------------------------
    %   gradiente real
    %----------------------------------------------------------------------
    grad = Btinv*grad;
    %----------------------------------------------------------------------
    %   calculo das constantes
    %----------------------------------------------------------------------
    wdt2jac = w(m)*dt2*jacob;
    %----------------------------------------------------------------------
    %   constantes da concentracao de poluicao
    %----------------------------------------------------------------------
    cc1 = wdt2jac*alfa(1);
    cc2 = wdt2jac*vx(1);
    cc3 = wdt2jac*vy(1);
    cc4 = dt2*sigma;
    %----------------------------------------------------------------------
    %   constantes da populacao p
    %----------------------------------------------------------------------
    cp1 = wdt2jac*alfa(2);
    cp2 = wdt2jac*vx(2);
    cp3 = wdt2jac*vy(2);    
    cp4 = dt2*lambda(2);
    %----------------------------------------------------------------------
    %   constantes da populacao q
    %----------------------------------------------------------------------
    cq1 = wdt2jac*alfa(3);
    cq2 = wdt2jac*vx(3);
    cq3 = wdt2jac*vy(3);
    cq4 = dt2*lambda(3);
    %----------------------------------------------------------------------
    for i = 1:3
        for j = 1:3
            %---=======================================--------------------
            %   submatrizes da concentracao de poluicao
            %---=======================================--------------------
            %   termo de difusao
            %--------------------------------------------------------------
            aux1 = cc1*(grad(:,j)')*grad(:,i);
            %--------------------------------------------------------------
            %   termo de adveccao
            %--------------------------------------------------------------
            aux2 = cc2*grad(1,j)*phi(i) + cc3*grad(2,j)*phi(i);
            %--------------------------------------------------------------
            %   termo de "reacao"
            %--------------------------------------------------------------
            aux3 = w(m)*phi(i)*phi(j)*jacob;
            %--------------------------------------------------------------
            %   soma de todos os termos
            %--------------------------------------------------------------
            meck(i,j) = meck(i,j) + aux1 + aux2 + (1 + cc4)*aux3;
            mdck(i,j) = mdck(i,j) - aux1 - aux2 + (1 - cc4)*aux3;
            %---==========================---------------------------------
            %   submatrizes da populacao p 
            %---==========================---------------------------------
            %   termo de difusao
            %--------------------------------------------------------------
            aux1 = cp1*(grad(:,j)')*grad(:,i);
            %--------------------------------------------------------------
            %   termo de adveccao
            %--------------------------------------------------------------
            aux2 = cp2*grad(1,j)*phi(i) + cp3*grad(2,j)*phi(i);
            %--------------------------------------------------------------      
            %   termo de "reacao" e' o mesmo
            %--------------------------------------------------------------
            %   soma de todos os termos
            %--------------------------------------------------------------
            mepk(i,j) = mepk(i,j) + aux1 + aux2 + (1 - cp4)*aux3;
            mdpk(i,j) = mdpk(i,j) - aux1 - aux2 + (1 + cp4)*aux3;
            %---==========================---------------------------------
            %   submatrizes da populacao q
            %---==========================---------------------------------
            %   termo de difusao
            %--------------------------------------------------------------
            aux1 = cq1*(grad(:,j)')*grad(:,i);
            %--------------------------------------------------------------
            %   termo de adveccao
            %--------------------------------------------------------------
            aux2 = cq2*grad(1,j)*phi(i) + cq3*grad(2,j)*phi(i);
            %--------------------------------------------------------------     
            %   termo de "reacao" e' o mesmo
            %--------------------------------------------------------------
            %   soma de todos os termos
            %--------------------------------------------------------------
            meqk(i,j) = meqk(i,j) + aux1 + aux2 + (1 - cq4)*aux3;
            mdqk(i,j) = mdqk(i,j) - aux1 - aux2 + (1 + cq4)*aux3;
            %--------------------------------------------------------------
        end
    end
end
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------