function B=try_BH(N,q)
% 2022-07-29
% random search for BH(N,q) for
% (15,12)
% (16,15)
% (17,12)
% (17,15)
% (18,14)
% (19,15)


    while(1);
        R=zeros(N);
	R(2:N,2:N)=(randi(q,N-1,N-1)-1);
	B=exp(2j*pi*R/q);
	n1=norm(B*B'-N*eye(N),'fro');
	if n1<1e-8, B, break, end;
    end

end


