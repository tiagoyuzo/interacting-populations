%--------------------------------------------------------------------------
%   script principal
%
%   utiliza o metodo de elementos finitos com base de lagrange linear
%   no espaco e crank-nicolson no tempo para resolver um sistema
%   de equacoes com uma populacao (com verhulst) e um poluente
%--------------------------------------------------------------------------
clear all
clc
tic
%---=====------------------------------------------------------------------
%   dados
%---=====------------------------------------------------------------------
%   variaveis globais
%--------------------------------------------------------------------------
global nx ny nnx nny nn L H nel nt dt
%--------------------------------------------------------------------------
%   numero do experimento para nome das figuras
%--------------------------------------------------------------------------
expnum = 40;
%--------------------------------------------------------------------------
%   numero de elementos em x e em y
%
%   nx: numero de divisoes em x
%   ny: numero de divisoes em y
%   nnx: numero de nos em x
%   nny: numero de nos em y
%   nn = numero total de nos
%--------------------------------------------------------------------------
nx = 16;
ny = nx;
nnx = nx + 1;
nny = ny + 1;
nn = nnx*nny;
%--------------------------------------------------------------------------
%   tamanho do passo no tempo e tempo final
%
%   dt: tamanho do passo no tempo
%   tf: tempo final
%   nt: numero de passos no tempo
%--------------------------------------------------------------------------
dt = 0.01;
tf = 1;
nt = tf/dt;
%--------------------------------------------------------------------------
%   tamanho do retangulo
%
%   L: largura do retangulo
%   H: altura do retangulo
%   dx: tamanho do passo na direcao x
%   dy: tamanho do passo na direcao y
%--------------------------------------------------------------------------
L = 0.5;
H = L;
dx = L/nx;
dy = H/ny;
%--------------------------------------------------------------------------
%   definicao das condicoes de contorno
%
%   tipo = 0 se eh neumann homogenea
%   tipo = 1 se eh dirichlet homogenea
%   tipo(i) = -i se eh mista
%
%               tipo(4)
%           ----------------
%          |                |
%   tipo(1)|                | tipo(3)
%          |                |
%           ----------------
%               tipo(2)
%--------------------------------------------------------------------------
tipo = zeros(1,4);
for i = 1:4
    tipo(i) = -i;
end
% tipo(1) = 1;
%--------------------------------------------------------------------------
%   construcao da malha
%--------------------------------------------------------------------------
[C,verx,very,inb,xt,yt,ind] = malha(tipo);
%--------------------------------------------------------------------------
%   matriz de conectividade
%--------------------------------------------------------------------------
[M,nel] = conect;
%---=============----------------------------------------------------------
%   fonte pontual
%---=============----------------------------------------------------------
%   f: vetor com a fonte pontual
%--------------------------------------------------------------------------
f = zeros(nn,1);
f(ind(5),1) = 2;
%---=========================----------------------------------------------
%   condicoes de estabilidade
%---=========================----------------------------------------------
[alfa,vx,vy,~,~,~,~,~,~,~] = parametros;
%--------------------------------------------------------------------------
%   numero de Peclet: estabilidade relacionando difusao e adveccao
%--------------------------------------------------------------------------
% for i = 1:3
%     if (abs(vx(i)*dx/alfa(i)) > 2 || abs(vy(i)*dy/alfa(i) > 2))
%         error('viola condicao de Peclet!')
%     end
% end
% %--------------------------------------------------------------------------
% %   estabilidade entre espaco e tempo
% %--------------------------------------------------------------------------
% for i = 1:3
%     if( alfa(i)*dt/(dx*dx) > 1)
%         error('viola condicao de estabilidade!')
%     end
% end
%---==================-----------------------------------------------------
%   condicoes iniciais
%---==================-----------------------------------------------------
[c,p,q,cst,pst,qst,caux,paux,qaux] = condini(xt);
%---============================-------------------------------------------
%   graficos da condicao inicial
%---============================-------------------------------------------
%   formato e armazenamento
%--------------------------------------------------------------------------
[verc,caux] = grafico(c,caux,xt,ind,0);
[verp,paux] = grafico(p,paux,xt,ind,0);
[verq,qaux] = grafico(q,qaux,xt,ind,0);
%--------------------------------------------------------------------------
%   construcao da figura
%--------------------------------------------------------------------------
grafsurf(verx,very,verc,verp,verq,expnum,0);
%--------------------------------------------------------------------------
pause(0.01)
%---=========================----------------------------------------------
%   sistema - parte constante
%---=========================----------------------------------------------
[mec,mdc,mep,mdp,meq,mdq,bc,bp,bq] = sistema(M,C,inb);
[Lc,Uc] = lu(mec);
%---===================----------------------------------------------------
%   iteracoes temporais
%   ===================
%   a cada iteracao:
%
%   resolve-se um sistema linear para c
%       o valor de c eh atualizado nele mesmo
%   resolve-se um sistema nao linear para p em 4 iteracoes internas
%       pst eh o valor de p^(n+1) a cada iteracao interna
%       ao fim das iteracoes internas p^n = pst
%   resolve-se um sistema nao linear para q como acima
%--------------------------------------------------------------------------
%   loop das iteracoes temporais
%--------------------------------------------------------------------------
for it = 1:nt
    %----------------==========================----------------------------
    %   resolucao de mec*c^(n+1) = mdc*c^n + bc
    %----------------==========================----------------------------
    A = mec;
    b = mdc*c + dt*f;
    %----------------------------------------------------------------------
    cst = Uc\(Lc\b);
    for jl = 1:nn
        if (cst(jl) < 10e-8)
            cst(jl) = 0;
        end
    end
    %----------------------------------------------------------------------
    %   loop das iteracoes intermediarias
    %----------------------------------------------------------------------
    for il = 1:4
        %---===========================------------------------------------
        %   inclusao da nao linearidade
        %---===========================------------------------------------
        %   construcao do sistema - parte nao linear
        %
        %   mepnl: matriz a esquerda de p nao linear
        %   mdpnl: matriz a direita de p nao linear
        %   bpnl: vetor do lado direito de p nao linear
        %   meqnl: matriz a esquerda de q nao linear
        %   mdqnl: matriz a direita de q nao linear
        %   bqnl: vetor do lado direito de q nao linear
        %------------------------------------------------------------------
        [mepnl,mdpnl,bpnl,meqnl,mdqnl,bqnl] = ...
            sistemanaolin(M,C,inb,c,cst,p,pst,q,qst);
        %----------------================================================--
        %   resolucao de (mep + mepnl)*p^(n+1) = (mdp + mdpnl)*p^n + bpnl
        %----------------================================================--
        A = mep + mepnl;
        b = (mdp + mdpnl)*p + bpnl;        
        %------------------------------------------------------------------
        pst = lu_solver(A,b);
        %----------------================================================--
        %   resolucao de (meq + meqnl)*p^(n+1) = (mdq + mdqnl)*q^n + bqnl
        %----------------================================================--
        A = meq + meqnl;
        b = (mdq + mdqnl)*q + bqnl;
        %------------------------------------------------------------------
        qst = lu_solver(A,b);
        %------------------------------------------------------------------
        for jl = 1:nn
            if (pst(jl) < 10e-8)
                pst(jl) = 0;
            end
            if (qst(jl) < 10e-8)
                qst(jl) = 0;
            end
        end
        %------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
    %   atualizacao da iteracao
    %----------------------------------------------------------------------
    c = cst;
    p = pst;
    q = qst;
    %----------------------------------------------------------------------
    %   armazenamento para os graficos temporais
    %----------------------------------------------------------------------
    [verc,caux] = grafico(c,caux,xt,ind,it);
    [verp,paux] = grafico(p,paux,xt,ind,it);
    [verq,qaux] = grafico(q,qaux,xt,ind,it);
    %---========-----------------------------------------------------------
    %   graficos - somente em algumas iteracoes
    %---========-----------------------------------------------------------
    if it == nt/4 || it == nt/2 || it == nt
        %------------------------------------------------------------------
        grafsurf(verx,very,verc,verp,verq,expnum,it);
        %------------------------------------------------------------------        
    end
    %----------------------------------------------------------------------
    %   acompanhamento das iteracoes
    %----------------------------------------------------------------------
    it
    pause(0.01)
end
%---===================----------------------------------------------------
%   graficos evolutivos
%---===================----------------------------------------------------
t = 0:dt:tf;
t = t';
%--------------------------------------------------------------------------
graftemp(t,caux,paux,qaux,it,nt,expnum);
%--------------------------------------------------------------------------
toc
%--------------------------------------------------------------------------
%   fim
%--------------------------------------------------------------------------
% save('p-sig4', 'p')
% save('q-sig4', 'q')
% save('c-sig4', 'c')