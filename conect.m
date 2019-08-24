%--------------------------------------------------------------------------
%   matriz de conectividade
%
%   M: matriz 3 x nel, onde nel eh o numero de elementos;
%   M(i,j) = i-esimo vertice do j-esimo elemento finito,
%            onde i = 1,...,3 e j = 1,...,nel
%
%   programa baseado no de Pulino, P., obtido em:
%   http://www.ime.unicamp.br/~pulino/MT623/programas/
%--------------------------------------------------------------------------
function [M,nel] = conect
%--------------------------------------------------------------------------
global nx ny nny
%--------------------------------------------------------------------------
nel = 0;
for j = 1:ny
    ind = (j - 1)*nny;
    for i = 1:nx
        k = nel + 1;
        
        M(1,k) = ind + i;
        M(2,k) = M(1,k) + 1;
        M(3,k) = M(2,k) + nx;
        
        nel = k + 1;
        
        M(1,nel) = M(3,k) + 1;
        M(2,nel) = M(3,k);
        M(3,nel) = M(2,k);    
    end
end
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------
