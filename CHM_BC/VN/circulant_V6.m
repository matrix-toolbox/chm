function [H, x, is_unimodular]=circulant_V6(x0)
% 2021-09-17
% usage:
% for j=1:1000; [H x u]=circulant_V9(exp(2j*pi*rand(1,9))); if abs(mean(u)-u(1))<0.0001;disp('FOUND'); break; end ,  end

        [x info]=fsolve(@uc, x0);
        is_unimodular=abs(x)
        is_convergent=sum(abs(info)) % should be all zeros
        K=[
 x(1) x(2) x(3) x(4) x(5) x(6) ;
 x(6) x(1) x(2) x(3) x(4) x(5) ;
 x(5) x(6) x(1) x(2) x(3) x(4) ;
 x(4) x(5) x(6) x(1) x(2) x(3) ;
 x(3) x(4) x(5) x(6) x(1) x(2) ;
 x(2) x(3) x(4) x(5) x(6) x(1) ;

        ]/abs(x(1));
	H=dephase(K);
	N=size(H)(1);
        spectrum=transpose(abs(eig(H)))(1)
        F_norm=norm(H*H'-N*eye(N),'fro')
        defect=ud(H,'S',0.000001)
        save('circulant_V6.dat', '*')
end

function y=uc(x)
 y(1)= x(1)/x(2) + x(2)/x(3) + x(3)/x(4) + x(4)/x(5) + x(5)/x(6) + x(6)/x(1) ;
 y(2)= x(1)/x(3) + x(2)/x(4) + x(3)/x(5) + x(4)/x(6) + x(5)/x(1) + x(6)/x(2) ;
 y(3)= x(1)/x(4) + x(2)/x(5) + x(3)/x(6) + x(4)/x(1) + x(5)/x(2) + x(6)/x(3) ;
 y(4)= x(1)/x(5) + x(2)/x(6) + x(3)/x(1) + x(4)/x(2) + x(5)/x(3) + x(6)/x(4) ;
 y(5)= x(1)/x(6) + x(2)/x(1) + x(3)/x(2) + x(4)/x(3) + x(5)/x(4) + x(6)/x(5) ;
end


