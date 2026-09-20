function H = GH_expand(file)
% ------------------------------------------------------------------------------
% 2026-09-18 Wojciech Bruzda, name[at]uj.edu.pl : name = w.bruzda, https://matrix-toolbox.github.io/chm
% ------------------------------------------------------------------------------
% Rebuild a bordered group-developed complex Hadamard matrix from its generator.
%
% A file GH_N_d_L_group.gen holds one row: the phases of F : G -> U(1), in the order of the group elements.
% The matrix is
%
%     H(*, *) = H(*, h) = H(g, *) = 1
%     H(g, h) = F(h g^-1),
%
% so the core is F read through the index table I(i,j) = position of g_j g_i^-1,
% which GH_groups.tsv stores for every group used here.
%
% >> H = GH_expand("GH_13_0_89_A4.gen");
% ------------------------------------------------------------------------------

    here = fileparts(mfilename("fullpath"));
    [~, stem] = fileparts(file);
    % GH_N_d_L_group, or GH_N_d_L_group__k when two classes share the whole
    % label.  Parsed with a regexp rather than strsplit, which collapses the
    % doubled underscore and would turn SD16__2 into the group name SD16_2 --
    % and Z8sZ4_2 really is a group name, so the two cannot be told apart
    % afterwards.
    tok = regexp(stem, '^GH_(\d+)_\d+_\d+_(.+?)(?:__\d+)?$', "tokens");
    if isempty(tok)
        error("GH_expand: '%s' is not a GH_N_d_L_group name", stem);
    end
    want = str2double(tok{1}{1});                % the N in the file name
    grp = tok{1}{2};

    a = load("-ascii", file);
    if any(~isfinite(a(:)))
        error("GH_expand: '%s' has non-finite phases", file);
    end
    F = exp(2i * pi * a(:)).';
    n = numel(F);
    if ~isfinite(want) || want ~= n + 1
        error("GH_expand: name says N=%g but the generator gives N=%d", want, n + 1);
    end

    I = [];
    fid = fopen(fullfile(here, "GH_groups.tsv"), "r");
    % assignment inside the condition is an Octave extension; written out so
    % the reader also runs under MATLAB
    line = fgetl(fid);
    while ischar(line)
        if ~isempty(line) && line(1) ~= "#"
            c = strsplit(line, "\t");
            if strcmp(c{1}, grp)
                I = reshape(sscanf(c{3}, "%d"), n, n).';
                break;
            end
        end
        line = fgetl(fid);
    end
    fclose(fid);
    if isempty(I)
        error("GH_expand: group '%s' not found in GH_groups.tsv", grp);
    end

    H = ones(n + 1);
    H(2:end, 2:end) = F(I);
end
