function [ pts ] = initialise_particles(flags, params, Np, observ)
%INITIALISE_PARTICLES Set up particles for particle filtering

ds = params.state_dim;
do = params.obs_dim;
dr = params.rnd_dim;

% Get position of first observation
if flags.space_dim == 2
    if flags.obs_mod == 1
        r_init = observ(1:2,1);
    elseif flags.obs_mod == 2
        [r1, r2] = pol2cart(observ(1,1), observ(2,1));
        r_init = [r1; r2];
    end
    
elseif flags.space_dim == 3
    if flags.obs_mod == 1
        r_init = observ(1:3,1);
    elseif flags.obs_mod == 2
        [r1, r2, r3] = sph2cart(observ(1,1), observ(2,1), observ(3,1));
        r_init = [r1; r2; r3];
    end
    
end

% Get velocity, or use a priori value if not observed
if flags.space_dim == 2
    if flags.obs_vel
        if flags.obs_mod == 1
            v_init = observ(3:4,1);
        elseif flags.obs_mod == 2
            [v1, v2] = pol2cart(observ(3,1), observ(4,1));
            v_init = [v1; v2];
        end
    else
        v_init = params.start_state(3:4);
    end
    
elseif flags.space_dim == 3
    if flags.obs_vel
        if flags.obs_mod == 1
            v_init = observ(4:6,1);
        elseif flags.obs_mod == 2
            [v1, v2, v3] = sph2cart(observ(4,1), observ(5,1), observ(6,1));
            v_init = [v1; v2; v3];
        end
    else
        v_init = params.start_state(4:6);
    end
    
end

x_init = [r_init; v_init];

% Generate a random set of starting points and accelerations
start_pos = mvnrnd(repmat(x_init', Np, 1), params.start_var)';
start_pos(4,:) = max(start_pos(4,:), params.min_speed);
start_accel = mvnrnd(zeros(Np,dr), params.Q)';
w_prob = log(mvnpdf(start_accel', zeros(Np,dr), params.Q));

% Put the values into cell array
pts = struct('x', num2cell(start_pos, [4,1])',...
             'w', num2cell(start_accel, [4,1])',...
             'tau', 0,...
             'Ns', 1,...
             'intx', num2cell(start_pos, [4,1])',...
             'lhood', 0, ...
             'tau_prob', 0,...
             'w_prob', num2cell(w_prob, 2));

end