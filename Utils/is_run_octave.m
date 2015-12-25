function [ result ] = is_run_octave( )
%IS_RUN_IN_MATLAB Check that the code running environment is Octave
%   See https://www.gnu.org/software/octave/doc/interpreter/How-to-distinguish-between-Octave-and-Matlab_003f.html

    persistent cache_result;  % speeds up repeated calls

    if isempty (cache_result)
        cache_result = exist('OCTAVE_VERSION', 'builtin') > 0;
    end

    result = cache_result;
end

