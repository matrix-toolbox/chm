function ER_pair(H)
% 2017-12-12
%
% Seek for ER rows (columns) of a CHM matrix.
%
% Call:
% >> er_pair(FN(4))
%
% Details:
%   D. Goyeneche, A new method to construct families of CHM in even dimensions
%   https://arxiv.org/abs/1210.7673
%   Def. III.1
%   Let CA and CB be two columns [rows] of a CHM. They are ER iff abs(real(CA*(j).CB(j)))=1 for any j.
%
%   Thm. III.1
%   Let H be a CHM defined in an even dimension d > 2. if H has m < d / 2 ER pairs,
%   then H belongs to a m-dim. family.

for r1=1:size(H, 1)-1
    for r2=r1+1:size(H, 1)
        ER = true;
        for j=1:size(H, 1)
            if abs(abs(real(H(r1, j)'*H(r2,j))) - 1) >= 1e-8
                ER = false;
            end
        end
        if ER == true
            printf("r1 = %d \t r2 = %d\n", r1, r2);
        end
    end
end

end

