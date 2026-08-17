function SL3(X, dSLMode)
% 20200411
% 20200927
% display a triplet of entropies NOT NORMALIZED TO UNITY!!

	try
		switch (dSLMode)
			case {"verbose"}
			disp(sprintf( 'SL(X,X^R,X^T) = %2.12g\t%2.12g\t%2.12g', SL(X), SL(R(X)), SL(T2(X)) ));
			case {"minimal"}
			
            if abs(SL(X)-SL(R(X)))<1e-10 && abs(SL(X) - SL(T(X)))<1e-10
                if abs(SL(X)-1.0)<1e-10
                    disp(sprintf( '%2.10g\t%2.10g\t%2.10g *', SL(X), SL(R(X)), SL(T2(X)) ));
                else
                    disp(sprintf( '%2.10g\t%2.10g\t%2.10g =', SL(X), SL(R(X)), SL(T2(X)) ));
                end
            else
                disp(sprintf( '%2.10g\t%2.10g\t%2.10g', SL(X), SL(R(X)), SL(T2(X)) ));
            end
            
			case {"compact"}
			disp(sprintf( '%2.12g\n%2.12g\n%2.12g', SL(X), SL(R(X)), SL(T2(X)) ));
			otherwise
			error('wrong display mode, use ''verbose'', ''minimal'' or ''compact''');
		end
	catch
		disp("missing parameter - default view:");
		disp(sprintf( '%2.12g\n%2.12g\n%2.12g', SL(X), SL(R(X)), SL(T2(X)) ));
	end

end


