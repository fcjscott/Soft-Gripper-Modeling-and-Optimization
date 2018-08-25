 
clear all;
clc;
%-----------------------------------------------

global miu object radius
object = "Cubic";
miu = 0.4;
angle = 70;

% initial guess of epsilon and delta
epsilon = 0.01;
delta = -0.002;
err1 = 1;
err2 = 1;
ct = 1;

while abs(err1) > 0.1 || abs(err2) > 0.01
    clf(figure(2));
    i = 1;
    % initial guess of finger spacing
    fingerSpacing = 0.018;
    verticalMove = -0.06;
    for radius = 0.02:0.005:0.025
        [Ps,Fy] = DER(fingerSpacing,verticalMove,angle);
        figure(2)
        plot(Ps(1:end-1),2*Fy,'o','MarkerSize',3);
        ylabel('Lifting Force (N)');
        xlabel('Instant Pressure (kPa)');
        title('Least square model for sphere','FontSize',14);
        grid on;
        hold on;
        legend({'radius = 0.02 m','radius = 0.025 m'},'FontSize',12);
        LiftArray(i,:) = 2*Fy;
        % update spacing and gripper position
        if i == 2
            break;
        end
        fingerSpacing = fingerSpacing + epsilon;
        verticalMove = verticalMove + delta;
        i = i + 1;
    end
    err1 = mean(LiftArray(1,:)-LiftArray(2,:));
    err2 = var(LiftArray(1,:)-LiftArray(2,:));
    disp(delta);
    disp(epsilon);
    pause(5);
    % update epsilon and delta based on error
    if err1 > 0.1 && err2 < 0.01
        epsilon = epsilon - 0.001;
    elseif err1 < -0.1 && err2 < 0.01
        epsilon = epsilon + 0.001;
    elseif err1 > 0.1 && err2 > 0.01
        epsilon = epsilon - 0.001;
        delta = delta - 0.001;
    elseif err1 < -0.1 && err2 > 0.01
        epsilon = epsilon - 0.001;
        delta = delta +0.001;
    end
    ct = ct + 1;
end

%%





