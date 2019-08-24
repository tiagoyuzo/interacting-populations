% %--------------------------------------------------------------------------
function [c,p,q,cst,pst,qst,caux,paux,qaux] = condini(xt)
%--------------------------------------------------------------------------
%   declaracao dos vetores com suas condicoes iniciais e vetores
%   auxiliares para construcao de graficos temporais
%
%   concentracao de poluicao - nula
%   populacao p - constante
%   populacao q - constante
%
%   pst e qst: populacoes auxiliares no tempo intermediario
%
%   caux, paux e qaux: poluicao auxiliar em nos fixos para graficos
%   temporais
%--------------------------------------------------------------------------
global nn nt
%--------------------------------------------------------------------------
%   c: concentracao de poluicao no tempo atual
%   caux: poluicao auxiliar em nos fixos para graficos temporais
%--------------------------------------------------------------------------
c = zeros(nn,1);
% load('c-sig8');
cst = c;
%--------------------------------------------------------------------------
% p = zeros(nn,1);
% p = 150*rand(nn,1);
p = 100*ones(nn,1);
% load('p-est');
% load('p-sig8');
pst = p;
%--------------------------------------------------------------------------
% q = zeros(nn,1);
% q = 150*rand(nn,1);
q = 100*ones(nn,1);
% load('q-est');
% load('q-sig8');
qst = q;
%--------------------------------------------------------------------------
qaux = zeros(nt+1,length(xt));
paux = zeros(nt+1,length(xt));
caux = zeros(nt+1,length(xt));
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------