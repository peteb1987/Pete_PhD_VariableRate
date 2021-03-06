function [ new_x ] = next_state( flags, params, old_x, w, dt )

% THIS IS THE OLD VERSION WHICH USES ROTATIONS.

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
    
    % Transform to planar intrinisics
    old_sdot = norm(old_v);
    
    if flags.space_dim == 2
        % Calculate 2D bearing
        old_psi = atan2(old_v(2), old_v(1));
        aNc = aN(1,:);
        
    elseif flags.space_dim == 3
        % Calculate declination and magnitude of normal acceleration (phi is
        % angle between the vector and a vertical plane.)
        [phi, aNc] = cart2pol(aN(1,:),aN(2,:));
        
        % Calculate unit vectors for 3D intrinsics at start of sojourn
        et = unit(old_v,1);
        et_rep = repmat(et,1,Ns);
        en = zeros(3,Ns); u = sqrt(et(1)^2+et(2)^2);
        en(1,:) =  (et(2)*sin(phi)-et(1)*et(3)*cos(phi))/u;
        en(2,:) = (-et(1)*sin(phi)-et(2)*et(3)*cos(phi))/u;
        en(3,:) = cos(phi)*u;
        eb = cross(et_rep,en);
        %     eb = zeros(3,Ns);
        %     part1 = (et(1)^2+et(2)^2);
        %     part2 = sqrt(cos(phi)^2-et(3)^2);
        %     eb(1,:) = (-et(1)*et(3)*sin(phi)+et(2)*part2)/part1;
        %     eb(2,:) = (-et(2)*et(3)*sin(phi)-et(1)*part2)/part1;
        %     eb(3,:) = sin(phi);
        %     en = cross(eb,et_rep);
        if Ns == 1
            R = [et, en, eb];
        else
            R = cat(3, et_rep', en', eb');
        end
        
        old_psi = 0;
        
    end
    
    %%% Solve planar differential equation %%%
    
    % speed
    aT = max(aT, (min_speed-old_sdot)/dt(end));
    new_sdot = old_sdot + aT*dt;
    
    % bearing
    if abs(aT)>1E-10
        new_psi = old_psi + (aNc./aT).*log(new_sdot./old_sdot);
    else
        new_psi = old_psi + (aNc.*dt)./old_sdot;
    end
    
    % displacement
    SF = 4*aT.^2 + aNc.^2;
    
    if Ns == 1
        
        new_u = zeros(flags.space_dim,no_cols);
        if (aT~=0)&&(aNc~=0)
            new_u(1,:) = ((new_sdot.^2)./SF).*( aNc.*sin(new_psi)+2*aT.*cos(new_psi)) - ((old_sdot^2)./SF)*( aNc.*sin(old_psi)+2*aT.*cos(old_psi));
            new_u(2,:) = ((new_sdot.^2)./SF).*(-aNc.*cos(new_psi)+2*aT.*sin(new_psi)) - ((old_sdot^2)./SF)*(-aNc.*cos(old_psi)+2*aT.*sin(old_psi));
        elseif (aT==0)&&(aNc~=0)
            new_u(1,:) = ((new_sdot.^2)./aNc).*( sin(new_psi) - sin(old_psi) );
            new_u(2,:) = ((new_sdot.^2)./aNc).*(-cos(new_psi) + cos(old_psi) );
        elseif (aT~=0)&&(aNc==0)
            new_u(1,:) = 0.5*dt.*cos(old_psi).*new_sdot;
            new_u(2,:) = 0.5*dt.*sin(old_psi).*new_sdot;
        else
            new_u(1,:) = ( old_sdot*dt.*cos(old_psi) );
            new_u(2,:) = ( old_sdot*dt.*sin(old_psi) );
        end
        
    elseif Ns>1
        
        new_u = zeros(flags.space_dim,no_cols);
        
        % Flags for where aT and aN are zero
        waTz = (aT==0);             % flag for where aT is 0
        waNz = (aNc==0);             % flag for where aN is 0
        bnz = ((~waTz)&(~waNz));   % both not zero
        oaTz = ((waTz)&(~waNz));   % only aT zero
        oaNz = ((~waTz)&(waNz));   % only aN zero
        bz = ((waTz)&(waNz));      % both zero
        
        % Neither aT nor aN are zero
        new_u(1,bnz) = ((new_sdot(bnz).^2)./SF(bnz)).*( aNc(bnz).*sin(new_psi(bnz))+2*aT(bnz).*cos(new_psi(bnz))) - ((old_sdot.^2)./SF(bnz)).*( aNc(bnz).*sin(old_psi)+2*aT(bnz).*cos(old_psi));
        new_u(2,bnz) = ((new_sdot(bnz).^2)./SF(bnz)).*(-aNc(bnz).*cos(new_psi(bnz))+2*aT(bnz).*sin(new_psi(bnz))) - ((old_sdot.^2)./SF(bnz)).*(-aNc(bnz).*cos(old_psi)+2*aT(bnz).*sin(old_psi));
        
        % aT is zero but aN isn't
        new_u(1,oaTz) = ((new_sdot(oaTz).^2)./aNc(oaTz)).*( sin(new_psi(oaTz)) - sin(old_psi) );
        new_u(2,oaTz) = ((new_sdot(oaTz).^2)./aNc(oaTz)).*(-cos(new_psi(oaTz)) + cos(old_psi) );
        
        % aN is zero but aT isn't
        new_u(1,oaNz) = 0.5*dt*cos(old_psi).*new_sdot(oaNz);
        new_u(2,oaNz) = 0.5*dt*sin(old_psi).*new_sdot(oaNz);
        
        % Both aT and aN are zero
        new_u(1,bz) = ( old_sdot*dt.*cos(old_psi) );
        new_u(2,bz) = ( old_sdot*dt.*sin(old_psi) );
        
    end
    
    % Calculate cartesian in-plane velocity
    new_udot = zeros(flags.space_dim,no_cols);
    [new_udot(1,:), new_udot(2,:)] = pol2cart(new_psi, new_sdot);
    
    % Translate back to original coordinate system
    if flags.space_dim == 2
        new_r = bsxfun(@plus, new_u+aX*dt, old_r);
        new_v = new_udot;
    elseif flags.space_dim == 3
        if Ns>1
            %         new_r = multiprod(R, new_u', [2,3], 2)' + aX*dt +repmat(old_r,1,Ns);
            %         new_v = multiprod(R, new_udot', [2,3], 2)';
            R_perm = permute(R,[2 3 1]);
            new_r = arraymatprod(R_perm, new_u) + aX*dt +repmat(old_r,1,Ns);
            new_v = arraymatprod(R_perm, new_udot);
        else
            new_r = bsxfun(@plus, R*new_u + aX*dt, old_r);
            new_v = R*new_udot;
        end
    end
    
end

% Stack them up
new_x = [new_r; new_v];

% Check
assert(all(isreal(new_x)));

end

