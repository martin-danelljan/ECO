function in_octave = is_octave()
%
% is_octave Test if run in GNU/Octave
%
% Return: true if the environment is Octave.
%
% Refer https://www.gnu.org/software/octave/doc/v4.0.1/How-to-distinguish-between-Octave-and-Matlab_003f.html
%

    persistent cacheval;  % speeds up repeated calls

    if isempty (cacheval)
        cacheval = (exist ('OCTAVE_VERSION', 'builtin') > 0);
    end

    in_octave = cacheval;
end