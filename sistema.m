%--------------------------------------------------------------------------
%   sistema - parte constante
%
%   constroi as submatrizes e subvetores, coloca as cond. de contorno
%   essenciais e constroi a parte linear do sistema
%
%   programa baseado no de Pulino, P., obtido em:
%   http://www.ime.unicamp.br/~pulino/MT623/programas/
%--------------------------------------------------------------------------
function [mec,mdc,mep,mdp,meq,mdq,bc,bp,bq] = sistema(M,C,inb)
%--------------------------------------------------------------------------
global nn nel
%--------------------------------------------------------------------------
%   declaracao usando a funcao sparse para economizar memoria
%   e operacoes
%
%   mec: matriz a esquerda de c
%   mdc: matriz a direita de c
%   mep: matriz a esquerda de p
%   mdp: matriz a direita de p
%   meq: matriz a esquerda de q
%   mdq: matriz a direita de q
%   bc: vetor do lado direito de c
%   bp: vetor do lado direito de p
%   bq: vetor do lado direito de q
%--------------------------------------------------------------------------
bc = sparse(nn,1);
bp = sparse(nn,1);
bq = sparse(nn,1);
mec = sparse(nn,nn);
mdc = sparse(nn,nn);
mep = sparse(nn,nn);
mdp = sparse(nn,nn);
meq = sparse(nn,nn);
mdq = sparse(nn,nn);
%--------------------------------------------------------------------------
%   loop de cada elemento: 1 ate nel
%--------------------------------------------------------------------------
for k = 1:nel
    %----------------------------------------------------------------------
    %   construcao da submatriz
    %----------------------------------------------------------------------
    [meck,mdck,mepk,mdpk,meqk,mdqk,bck,bpk,bqk] = submatriz(M,C,k);
    %----------------------------------------------------------------------
    %   condicao de contorno (essencial)
    %----------------------------------------------------------------------
%     [meck,mdck,bck] = cond_contorno(meck,mdck,bck,C,M,k,inb,1,1);
%     [mepk,mdpk,bpk] = cond_contorno(mepk,mdpk,bpk,C,M,k,inb,2,1);
%     [meqk,mdqk,bqk] = cond_contorno(meqk,mdqk,bqk,C,M,k,inb,3,1);
    %----------------------------------------------------------------------
    %   construcao do sistema
    %
    %   utiliza a funcao sparse e a matriz de conectividade
    %   cada elemento da submatriz eh colocado na matriz por meio da
    %   matriz de conectividade, elemento a elemento
    %----------------------------------------------------------------------
    for i = 1:3
        bc = bc + sparse(M(i,k),1,bck(i,1),nn,1);
        bp = bp + sparse(M(i,k),1,bpk(i,1),nn,1);
        bq = bq + sparse(M(i,k),1,bqk(i,1),nn,1);
        for j = 1:3
            mec = mec + sparse(M(i,k),M(j,k),meck(i,j),nn,nn);
            mdc = mdc + sparse(M(i,k),M(j,k),mdck(i,j),nn,nn);
            mep = mep + sparse(M(i,k),M(j,k),mepk(i,j),nn,nn);
            mdp = mdp + sparse(M(i,k),M(j,k),mdpk(i,j),nn,nn);
            meq = meq + sparse(M(i,k),M(j,k),meqk(i,j),nn,nn);
            mdq = mdq + sparse(M(i,k),M(j,k),mdqk(i,j),nn,nn);
        end
    end
end
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------