function get_defect
# 2022-07-30

    clear all;
    B=OCTAVE_INPUT;
    B_SIZE=size(B, 2);
    for k=1:B_SIZE;
        disp(sprintf("%g", ud(B{k}, "S", 1e-8)));
    end

end

