%--------------------------------------------------------------------------
%
%   verx e very: vetores que contem os valores de x e y respectivamente
%   malha
%
%   constroi a malha de elementos finitos de primeira ordem para um 
%   dominio retangular de primeira ordem, faz a separacao dos bordos
%   segundo diferentes condicoes de contorno e cria pontos auxiliares
%   para construcao dos graficos temporais
%
%   a enumeracao dos nos e' feita da esquerda para direita
%
%   programa baseado no de Pulino, P., obtido em:
%   http://www.ime.unicamp.br/~pulino/MT623/programas/
%--------------------------------------------------------------------------
function [C,verx,very,inb,xt,yt,ind] = malha(tipo)
%--------------------------------------------------------------------------
global nx ny nnx nny nn L H
%--------------------------------------------------------------------------
%   tamanho do passo
%
%   dx: tamanho do passo em x
%   dt: tamanho do passo em y
%--------------------------------------------------------------------------
dx = L/nx;
dy = H/ny;
%--------------------------------------------------------------------------
%   vetores de saida de x e y
%--------------------------------------------------------------------------
verx = zeros(nnx,1);
very = zeros(nny,1);

for i = 2:nnx
    verx(i) = verx(i-1) + dx;
end

for i = 2:nny
    very(i) = very(i-1) + dy;
end
%--------------------------------------------------------------------------
%   matriz de coordenadas
%
%   C: matriz 2 x nn com as coordenadas de x na primeira linha e de y
%      na segunda; cada coluna representa um no
%--------------------------------------------------------------------------
k = 0;
for j = 1:nny
    for i = 1:nnx
        k = k + 1;
        C(1,k) = verx(i);
        C(2,k) = very(j);
    end
end
%--------------------------------------------------------------------------
%   identificacao dos nos do contorno
%
%   inb > 0 se a c.c. eh dirichlet
%   inb < 0 se a c.c. eh mista
%   inb = 0 se a c.c. eh neumann homogenea
%
%               tipo(4)
%           ----------------
%          |                |
%   tipo(1)|                | tipo(3)
%          |                |
%           ----------------
%               tipo(2)
%--------------------------------------------------------------------------
inb = zeros(k,1);

for i = 1:nnx
    ind = i + nnx*ny;
    
    inb(i) = tipo(2);
    inb(ind) = tipo(4);
end

for i = 1:(ny-1)
    ind1 = i*nnx + 1;
    ind2 = ind1 + nx;
    
    inb(ind1) = tipo(1);
    inb(ind2) = tipo(3);
end
%--------------------------------------------------------------------
%   pontos para graficos temporais
% 
%   posicao dos pontos na malha
%
%    ----------------
%   |  1    2    3   |
%   |  4    5    6   |
%   |  7    8    9   |
%    ----------------
%--------------------------------------------------------------------
nosx = L*[0.25;0.5;0.75];
nosy = H*[0.25;0.5;0.75];
ind = zeros(9,1);
for i = 1:nn
    for j = 1:3
        for k = 1:3
            if((abs(C(1,i)-nosx(j,1)) < 10^-4) && ...
                    (abs(C(2,i)-nosy(k,1)) < 10^-4))
                ind(3*(k-1)+j) = i;
            end
        end
    end
end
xt = zeros(length(ind),1);
yt = zeros(length(ind),1);
for i = 1:length(ind)
    xt(i,1) = C(1,ind(i));
    yt(i,1) = C(2,ind(i));
end
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------