function timeIntegration(ki0)

global q q0
global nv dt
global mapCons ForceAll error
global uProjected

%% Step 1: Predictor
fprintf('Predictor step\n');
[q, error] = objfunBW(q0,ki0);
if (error < 0)
    fprintf('Could not converge. Sorry\n');
    return;
end

%% Detect and delete constrained dof
changeMade = 0; % flag to identify if another step is needed


for c=3:nv
    % Detect new constraints
    xPos = q(2*c-1);
    yPos = q(2*c);
    boundaryY = getConstraint(xPos,yPos);
    
    if ( mapCons(2*c) == 0 && abs(yPos+0.15) < boundaryY )
        
        mapCons(2*c) = 1; % constrain it

        % Compute projected velocity. This will be used in objfunBW to
        % impose constraints
        uProjected(2*c-1) = (q(2*c-1) - q0(2*c-1))/dt;
        uProjected(2*c) = (q(2*c) - q0(2*c))/dt;
        
        changeMade = 1;
    end
end

%% Step 2: Corrector for new constrained dofs
if (changeMade == 1)
    fprintf('Corrector-1 step\n');
    [q, error] = objfunBW(q0,ki0);
    if (error < 0)
        fprintf('Could not converge. Sorry\n');
        return;
    end    
end

%% Detect and delete constrained dof
changeMade = 0; % flag to identify if another step is needed

for c=3:nv
    % Detect new constraints
    xPos = q(2*c-1);
    yPos = q(2*c);
    boundaryY = getConstraint(xPos,yPos);

    % Delete unnecessary constraints
    if (mapCons(2*c) == 1)
        
        % Compute reaction force from ground
        fReaction = ForceAll(2*c-1:2*c);
        fNormal = dot(fReaction, nConstraint(xPos,yPos));
        
        % Based on reaction force, release the constraint
        if (fNormal <= 0) % reaction force is negative
            mapCons(2*c) = 0; % Unconstrain it
            changeMade = 1;
        end
    end
    
    % Detect new constraints
    if ( abs(yPos+0.15) < boundaryY )
        
        mapCons(2*c) = 1; % constrain it

        % Compute projected velocity. This will be used in objfunBW to
        % impose constraints
        uProjected(2*c-1) = (q(2*c-1) - q0(2*c-1))/dt;
        uProjected(2*c) = (q(2*c) - q0(2*c))/dt;
        
        changeMade = 1;
    end
end

%% Step 3: Corrector for released dofs
if (changeMade == 1)
    fprintf('Corrector-1 step\n');
    [q, error] = objfunBW(q0,ki0);
    if (error < 0)
        fprintf('Could not converge. Sorry\n');
        return;
    end    
end

end

