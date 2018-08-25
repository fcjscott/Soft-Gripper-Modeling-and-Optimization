function rod_id(x, ctime, saveImage, imageDirectory)
global dt object radius
x1 = x(1:2:end,1);
x2 = x(2:2:end,1);


h1 = figure(1);
%set(h1, 'visible', 'off');
clf()
plot(x1,x2, 'ko-');
hold on
plot(-x1,x2, 'ko-');
hold on
%Mark the clamped nodes red
plot(x1(1),x2(1), 'r^');
plot(x1(2),x2(2), 'r^');
plot(-x1(1),x2(1), 'r^');
plot(-x1(2),x2(2), 'r^');
hold on


%plot the object
if object == "Ball"
    theta = 0:pi/180:2*pi;
    r = radius;
    a = r*cos(theta);
    b = r*sin(theta)-0.15;
    plot(a,b);
elseif object == "Cubic"
    a1 = -radius;
    a2 = radius;
    b1 = -0.15+radius;
    b2 = -0.15-radius;
    a = [a1, a2, a2, a1, a1];
    b = [b1, b1, b2, b2, b1];
    plot(a,b,'-');
else
end


hold off
title(num2str(ctime, 't=%f \n'));
axis equal
xlabel('x [m]');
ylabel('y [m]');
grid on
drawnow
if (saveImage~=0)
    saveas(h1, num2str(ctime/dt, [imageDirectory, '/t=%09.0f.png']));
end

end