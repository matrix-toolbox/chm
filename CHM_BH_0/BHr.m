#!/usr/bin/octave -qf

arg_list = argv()

chunk = str2num(arg_list{1});
chunkNo = str2num(arg_list{2});
m =  str2num(arg_list{3});
q =  str2num(arg_list{4});
getButson = str2func(arg_list{5});

for j = 1:chunk
    Butson = getButson(j);
    s = qSymmetry2(Butson, m, q);
    if s<1e-10
%        BNumber = chunkNo*chunk + j;
        BNumber = j;
        d=ud(Butson);
        fileName = strcat('BH', int2str(m*q), '_', int2str(q), '_', int2str(d), '_', int2str(BNumber), '.dat')
        LButson = round(log(Butson)*q/2j/pi);
        save(fileName, 'LButson');
    end
end




