function err = rms_error(wave,waveref,kwargs)
%RMS_ERROR Calculates the rmse between two waves, assuming the same discretisation in x-direction
% INPUTS:
%   wave: solutionObject
%   waveref: solutionObject, reference wave to calculate the rmse with
%   optional t: double, if specified then both waves will be reduced to times t
%      (by default, the time array will be reduced to that of waveref)
% OUTPUTS:
%   none
arguments
    wave solutionObject
    waveref solutionObject
    kwargs.t double = []
end

if isempty(kwargs.t)
    kwargs.t = waveref.t;
end

wave = wave.reduceToTimes(kwargs.t);
waveref = waveref.reduceToTimes(kwargs.t);

err = rmse(wave.u, waveref.u, 1);

end