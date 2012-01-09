function [ pos_rmse, vel_rmse ] = check_smooth_error( flags, params, true_intx, smooth_pts )
%CHECK_SMOOTH_ERROR Calculate mean smoothed trajectory and error from true
%values

sd = flags.space_dim;

% Grab the particle interpolated states
pts_intx = cat(3, smooth_pts.intx);

% Average them
smooth_intx = mean(pts_intx, 3);

% Calculated position and velovity errors
pos_error = sqrt(sum( (true_intx(1:sd,:)-smooth_intx(1:sd,:)).^2, 1));
vel_error = sqrt(sum( (true_intx(sd+1:2*sd,:)-smooth_intx(sd+1:2*sd,:)).^2, 1));

% Calculate RMSEs
pos_rmse = sqrt(mean(pos_error.^2));
vel_rmse = sqrt(mean(vel_error.^2));

figure, plot(pos_error);
figure, plot(vel_error);



end

