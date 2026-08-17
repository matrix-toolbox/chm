function K9_D_circles()
% 2023-02-16
%
% Draw four circles inside the area of D (for K9_2).
% They represent for inner orbits of K9_2 with two different defects: 2 and 10.
%
% >> K9_D_circles;

    hold all;

    t=0:0.001:1;scatter(2*cos(2*pi*t), sqrt(3)+2*sin(2*pi*t), 1, "k", "filled")  % CU
    t=0:0.001:1;scatter(2*cos(2*pi*t), -sqrt(3)+2*sin(2*pi*t), 1, "k", "filled") % CD
    t=0:0.001:1;scatter(1+2*cos(2*pi*t), 2*sin(2*pi*t), 1, "k", "filled")        % CR
    t=0:0.001:1;scatter(-1+2*cos(2*pi*t), 2*sin(2*pi*t), 1, "k", "filled")       % CL

end


