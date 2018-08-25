function [q, error] = objfunBW(q,ki0)
%%
global dt q0 u m ForceAll nv Len
global ScaleSolver tol maximum_iter
global mMatInv Ident2nv imposedAcceleration mapCons
global uProjected object radius

%% Figure out the imposed acceleration
for c=1:nv
    xPos = q0(2*c-1);
    yPos = q0(2*c);
    % Case 1: both dofs are constrained
    if (mapCons(2*c-1)==1 && mapCons(2*c)==1)
        
        % Compute the velocity perpendicular to the ground
        uProjectedPoint = uProjected(2*c-1:2*c);
        uPerp = dot(nConstraint(xPos,yPos), uProjectedPoint) *...
        nConstraint(xPos,yPos);
        normU = norm(uPerp);
        if (normU < (Len/9.81)*1e-3 ) % very small
            normUinv = 0;
        else
            normUinv = 1 / normU;
        end
        
        if object == "Ball"
            Rq = sqrt(q0(2*c-1)^2+(q0(2*c)+0.15)^2);
            dRDesired = (Rq - radius);
        elseif object == "Cubic"
            Rq = q0(2*c-1);
            dRDesired = -Rq - radius;
        end
        tpseudo = dRDesired * normUinv;
        
        q0Point = [q0(2*c-1); q0(2*c)];
        u0Point = [u(2*c-1); u(2*c)];
        qDesired = q0Point + tpseudo*uProjectedPoint;
        
        imposedAcceleration(2*c-1:2*c) = (qDesired - q0Point)/dt - u0Point;
        mMatInv(2*c-1:2*c,2*c-1:2*c) = invMass(xPos,yPos) * [0 0; 0 0];        

    % Case 2: one dof is constrained
    elseif (mapCons(2*c)==1)
        if object == "Ball"
            Rq = sqrt(q0(2*c-1)^2+(q0(2*c)+0.15)^2);
            dRDesired = (Rq - radius);
        elseif object == "Cubic"
            Rq = q0(2*c-1);
            dRDesired = -Rq - radius;
        end
        q0Point = [q0(2*c-1); q0(2*c)];
        u0Point = [u(2*c-1); u(2*c)];
        qDesired = q0Point - nConstraint(xPos,yPos) * dRDesired;
        
        imposedAcceleration(2*c-1:2*c) = (qDesired - q0Point)/dt - ...
            nConstraint(xPos,yPos) * dot(u0Point, nConstraint(xPos,yPos));
        mMatInv(2*c-1:2*c,2*c-1:2*c) = invMass(xPos,yPos) * [1/m(2*c-1) 0; 0 1/m(2*c)];
    else
        imposedAcceleration(2*c-1:2*c) = 0;
        mMatInv(2*c-1:2*c,2*c-1:2*c) = [1/m(2*c-1) 0; 0 1/m(2*c)];
    end
end

%%
% Newton-Raphson scheme
iter = 0; % number of iterations
normf = tol*ScaleSolver*10; % norm of function value (initialized to a value higher 
% than tolerance) 
error = 1; % Start with a 'good' simulation (error=1 means no error)

% Initial guess for delta V
dV = (q - q0)/dt - u;

while (normf > tol*ScaleSolver)
    % Figure out the velocities
    uNew = u + dV;
    
    % Figure out the positions    
    q = q0 + dt * uNew;
    
    % Get forces
    [Fb, Jb] = getFb(q,ki0);
    [Fs, Js] = getFs(q);
    Fg = getFg(q); % external force
    
    Forces = (Fb + Fs + Fg);
    
    % Equation of motion
    ForceAll = m .* dV/dt - Forces; % actual force
    f = dV - dt * mMatInv * Forces - imposedAcceleration; % force used for Baraff-Witkin mass modification
    
   
    % Get the norm
    normf = norm( f ) * mean(m)/dt;

    if (normf > tol*ScaleSolver)
        % Manipulate the Jacobians
        Jelastic = Jb + Js;
        
        % Newton's update
        J = Ident2nv - dt^2 * mMatInv * Jelastic;
        
        dV = dV - J \ f;
        
        % Update iteration number
        iter = iter+1;
        
        if (iter>maximum_iter)
            error = -1; % return with an error signal
            return
        end
    end
end

fprintf('Total iterations=%d, error=%e\n', iter, normf);

end
