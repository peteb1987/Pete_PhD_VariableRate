function [ mu ] = tracking_calc_obs_mean( flags, x )
%TRACKING_CALC_OBS_MEAN Calculate the Gaussian mean of the observation
%distribution for various observation models

if flags.obs_mod == 1
    mu = x(1:2);
elseif flags.obs_mod == 1
    error('You''ve not written this yet');
end

end

