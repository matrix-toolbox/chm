function K9_D_hull()
% 2023-02-16
%
% Draw a hull of D for for K9_2.
%
% >> K9_D_hull;

    hold all;
    axis([-4 4 -4 4], "square");

    t=-atan(sqrt(15))/2/pi:0.0005:atan(sqrt(15))/2/pi;scatter(-1+4*cos(2*pi*t), 4*sin(2*pi*t), 1, "k", "filled")
    t=1/2-atan(sqrt(15))/2/pi:0.0005:1/2+atan(sqrt(15))/2/pi;scatter(+1+4*cos(2*pi*t), 4*sin(2*pi*t), 1, "k", "filled")
    


end


