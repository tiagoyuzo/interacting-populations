%--------------------------------------------------------------------------
%   sistema nao linear
%
%   calcula os termos nao lineares a serem colocados no sistema
%--------------------------------------------------------------------------
function [mepnl,mdpnl,bpnl,meqnl,mdqnl,bqnl] = ...
    sistemanaolin(M,C,inb,c,cst,p,pst,q,qst)
%--------------------------------------------------------------------------
global nn nel
%--------------------------------------------------------------------------
%   declaracao usando a funcao sparse para economizar memoria
%   e operacoes
%
%   mepnl: matriz a esquerda de p nao linear
%   mdpnl: matriz a direita de p nao linear
%   meqnl: matriz a esquerda de q nao linear
%   mdqnl: matriz a direita de q nao linear
%   bpnl: vetor do lado direito de p nao linear
%   bqnl: vetor do lado direito de q nao linear
%--------------------------------------------------------------------------
mepnl = sparse(nn,nn);
mdpnl = sparse(nn,nn);
meqnl = sparse(nn,nn);
mdqnl = sparse(nn,nn);
bpnl = sparse(nn,1);
bqnl = sparse(nn,1);
bpk = zeros(3,1);
bqk = zeros(3,1);
%--------------------------------------------------------------------------
%   loop de cada elemento: 1 ate nel
%--------------------------------------------------------------------------
for k = 1:nel
    %----------------------------------------------------------------------
    %   construcao da submatriz
    %----------------------------------------------------------------------
    [mepk,mdpk,meqk,mdqk] = submatnaolin(M,C,k,c,cst,p,pst,q,qst);
    %----------------------------------------------------------------------
    %   condicao de contorno
    %----------------------------------------------------------------------
%     [mepk,mdpk,bpk] = cond_contorno(mepk,mdpk,bpk,C,M,k,inb,2,0);
%     [meqk,mdqk,bqk] = cond_contorno(meqk,mdqk,bqk,C,M,k,inb,3,0);
%     [mepk,mdpk,~] = cond_contorno(mepk,mdpk,bpk,C,M,k,inb,2,0);
%     [meqk,mdqk,~] = cond_contorno(meqk,mdqk,bqk,C,M,k,inb,3,0);
    %----------------------------------------------------------------------
    %   construcao do sistema
    %
    %   utiliza a funcao sparse e a matriz de conectividade
    %   cada elemento da submatriz eh colocado na matriz por meio da
    %   matriz de conectividade, elemento a elemento
    %----------------------------------------------------------------------
    for i = 1:3
%         bpnl = bpnl + sparse(M(i,k),1,bpk(i,1),nn,1);
%         bqnl = bqnl + sparse(M(i,k),1,bqk(i,1),nn,1);
        for j = 1:3
            mepnl = mepnl + sparse(M(i,k),M(j,k),mepk(i,j),nn,nn);
            mdpnl = mdpnl + sparse(M(i,k),M(j,k),mdpk(i,j),nn,nn);
            meqnl = meqnl + sparse(M(i,k),M(j,k),meqk(i,j),nn,nn);
            mdqnl = mdqnl + sparse(M(i,k),M(j,k),mdqk(i,j),nn,nn);
        end
    end
end
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------