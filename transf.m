%--------------------------------------------------------------------------
%   transf
%
%   faz a transformacao afim do triangulo padrao para o triangulo
%   local, constroi a matriz jacobiana B, (B^t)^-1 e x, y, coordenadas
%   reais
%
%   programa baseado no de Pulino, P., obtido em:
%   http://www.ime.unicamp.br/~pulino/MT623/programas/
%--------------------------------------------------------------------------
function [Btinv,jacob] = transf(M,C,k)
%--------------------------------------------------------------------------
%   identificacao das coordenadas
%--------------------------------------------------------------------------
ip = M(1,k);
jp = M(2,k);
kp = M(3,k);
%--------------------------------------------------------------------------
%   matriz jacobiana
%--------------------------------------------------------------------------
B(1,:) = [( C(1,jp) - C(1,ip) )  ( C(1,kp) - C(1,ip) )];
B(2,:) = [( C(2,jp) - C(2,ip) )  ( C(2,kp) - C(2,ip) )];
%--------------------------------------------------------------------------
%   determinante
%--------------------------------------------------------------------------
jacob = B(1,1)*B(2,2) - B(1,2)*B(2,1);
%--------------------------------------------------------------------------
%   inversa transposta do jacobiano
%--------------------------------------------------------------------------
% Btinv = inv(B');
Btinv = 1/jacob*[B(2,2), -B(2,1); -B(1,2), B(1,1)];
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------