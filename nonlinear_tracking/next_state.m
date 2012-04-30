function [ new_x ] = next_state( flags, params, old_x, w, dt )
%NEXT_STATE Deterministically calculate the next state given the previous
%state, the random variable, and the time difference using the various
%dynamic models

% Two types of batch operation can be handled: time batch, for which dt is
% a row vector, and acceleration batch, for which w is a matrix of
% accelerations in columns.

K = numel(dt);
Ns = size(w,2);
assert((Ns==1)||(K==1), 'Cannot handle batches in both time and acceleration');
no_cols = max(Ns, K);

% If no time has passed (i.e. we're looking at an observation exactly
% after a jump) then return straight away. (This happens at t=0)
if (K==1)&&(dt==0)
    new_x=old_x;
    return
end

if flags.dyn_mod == 3
    
    % Cartesian case
    
    % Get accelerations
    a = w(1:flags.space_dim,:);
    
    % Get old state
    old_r = old_x(1:flags.space_dim,:);
    old_v = old_x(flags.space_dim+1:end,:);
    
    % Solve equations
    new_r = bsxfun(@plus, bsxfun(@plus, old_r, old_v * dt), 0.5 * a * dt.^2);
    new_v = bsxfun(@plus, old_v, a * dt);
    
else
    
    % Intrinsic case
    
    % Set minimum speed
    min_speed = params.min_speed;
    
    % Get accelerations
    aT = w(1,:);
    aN = w(2:flags.space_dim,:);
    if flags.dyn_mod == 2
        aX = w(flags.space_dim+1:end,:);
    else
        aX = zeros(flags.space_dim,Ns);
    end
    
    % Get old state
    old_r = old_x(1:flags.space_dim,:);
    old_v = old_x(flags.space_dim+1:end,:);
    
    % Calculate declination and magnitude of initial normal acceleration 
    % (phi is angle between the vector and a vertical plane.)
    if flags.space_dim == 3
        [phi, aNc] = cart2pol(aN(1,:),aN(2,:));
    elseif flags.space_dim == 2
        aNc = aN; phi = pi/2;
    end
    
    % Flags for where aT and aN are zero
    tol = 1E-10;
    waTz = (abs(aT)<tol);      % flag for where aT is 0
    waNz = (abs(aNc)<tol);     % flag for where aN is 0
    bnz = ((~waTz)&(~waNz));   % both not zero
    oaTz = ((waTz)&(~waNz));   % only aT zero
    oaNz = ((~waTz)&(waNz));   % only aN zero
    bz = ((waTz)&(waNz));      % both zero
    
    % Inital speed
    old_sdot = norm(old_v);
    aT = max(aT, (min_speed-old_sdot)/dt(end));
    
    % Unit vector matrixes
    
    % Tangential
    et = unit(old_v,1);
    et_rep = repmat(et,1,Ns);
    
    % Normal
    if flags.space_dim == 3
        en = zeros(3,Ns); u = sqrt(et(1)^2+et(2)^2);
        en(1,:) =  (et(2)*sin(phi)-et(1)*et(3)*cos(phi))/u;
        en(2,:) = (-et(1)*sin(phi)-et(2)*et(3)*cos(phi))/u;
        en(3,:) = cos(phi)*u;
    elseif flags.space_dim == 2
        en = zeros(3,Ns);
        en(1,:) = et(2);
        en(2,:) = -et(1);
    end
    
    % Binormal
    eb = cross(et_rep,en);
    
    % Transitions
    
    % Speed
    new_sdot = old_sdot + aT*dt;
    
    % Angle
    dpsi = zeros(1, no_cols);
    if any(bnz)
        dpsi(bnz&true(1,K)) = (aNc(bnz)./aT(bnz)).*log(new_sdot(bnz&true(1,K))./old_sdot);
    end
    if any(oaTz)
        dpsi(oaTz&true(1,K)) = (aNc(oaTz).*dt)./old_sdot;
    end
    
    % Unit vectors
    new_et =  bsxfun(@times, cos(dpsi), et_rep) + bsxfun(@times, sin(dpsi), en);
    new_en = -bsxfun(@times, sin(dpsi), et_rep) + bsxfun(@times, cos(dpsi), en);
    new_eb = eb;
    
    % Convert back to cartesians
    
    % Velocity
    new_v = bsxfun(@times, new_sdot, new_et);
    
    % Displacement
    new_r = zeros(size(new_v));
    interm1 = bsxfun(@times, new_sdot(~bz&true(1,K)).^2, 2*bsxfun(@times, aT(~bz), cos(dpsi(~bz&true(1,K)))) + bsxfun(@times, aNc(~bz), sin(dpsi(~bz&true(1,K))))) - 2*aT(~bz)*old_sdot^2;
    interm2 = bsxfun(@times, new_sdot(~bz&true(1,K)).^2, 2*bsxfun(@times, aT(~bz), sin(dpsi(~bz&true(1,K)))) - bsxfun(@times, aNc(~bz), cos(dpsi(~bz&true(1,K))))) +  aNc(~bz)*old_sdot^2;
    
    new_r(:,~bz&true(1,K)) = bsxfun(@plus, bsxfun(@plus, old_r, aX(:,~bz)*dt), bsxfun(@times, ( 1 ./ ( aNc(~bz).^2+4*aT(~bz).^2 ) ), ...
        bsxfun(@times, interm1, et_rep(:,~bz)) + ...
        bsxfun(@times, interm2, en(:,~bz)) ) );
    if any(bz)
        new_r(:,bz&true(1,K)) = bsxfun(@plus, bsxfun(@plus, old_r, aX(:,bz)*dt), bsxfun(@times, new_v(:,bz&true(1,K)), dt) );
    end
    
end
    
% Stack them up
new_x = [new_r; new_v];

% Check
assert(all(isreal(new_x)));

end

