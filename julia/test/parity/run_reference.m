% run_reference.m
%
% Load inputs_<fn>.mat (flat schema: ncases, nargs, case<i>_arg<k>) and call
% the ORIGINAL SPM .m function from the repo root, writing each output to
% reference_<fn>.mat as out_<i> column vectors.
%
%   octave-cli --eval "run_reference('Npdf')"

function run_reference(fn)
    here = fileparts(mfilename('fullpath'));
    spm_root = fullfile(here, '..', '..', '..');
    addpath(spm_root);
    data_dir = fullfile(here, 'data');

    in_file  = fullfile(data_dir, sprintf('inputs_%s.mat', fn));
    out_file = fullfile(data_dir, sprintf('reference_%s.mat', fn));

    S = load(in_file);
    ncases = S.ncases;
    nargs  = S.nargs;

    warning('off', 'all');
    fhandle = str2func(sprintf('spm_%s', fn));

    out_struct = struct();
    out_struct.ncases = ncases;
    for i = 1:ncases
        args = cell(1, nargs);
        for k = 1:nargs
            args{k} = S.(sprintf('case%d_arg%d', i, k));
        end
        y = fhandle(args{:});
        out_struct.(sprintf('out_%d', i)) = y(:);   % flatten column
    end

    save('-v7', out_file, '-struct', 'out_struct');
    fprintf('wrote %s with %d outputs\n', out_file, ncases);
end
