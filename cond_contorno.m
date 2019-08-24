%--------------------------------------------------------------------------
%   condicao de contorno
%
%   aplica as condicoes de contorno em cada submatriz do sistema
%
%   programa baseado no de Pulino, P., obtido em:
%   http://www.ime.unicamp.br/~pulino/MT623/programas/
%--------------------------------------------------------------------------
function [mek,mdk,bk] = cond_contorno(mek,mdk,bk,C,M,k,inb,pop,linear)
%--------------------------------------------------------------------------
im = zeros(3,1);
for i = 1:3
    im(i,1) = M(i,k);
end
%--------------------------------------------------------------------------
for i = 1:3
    %----------------------------------------------------------------------
    %   se a condicao e' dirichlet homogenea
    %----------------------------------------------------------------------
    if( inb( im(i) ) > 0 )
        bk(i) = 0;
        for k = 1:3
            mek(k,i) = 0;
            mek(i,k) = 0;
            mdk(k,i) = 0;
            mdk(i,k) = 0;
        end
        mek(i,i) = 1;
        %------------------------------------------------------------------
        %   se a condicao e' de robin
        %------------------------------------------------------------------
    elseif(inb(im(i)) < 0 && linear == 1)
        for j = 1:3
            if(inb(im(j)) < 0 && i ~= j)
%             if(inb(im(j)) == inb(im(i)))
                [mek,mdk] = mixcont(inb(im(i)),im(i),im(j),mek,mdk,C,pop);
%                 return
            end
        end
        return
    end
end
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------