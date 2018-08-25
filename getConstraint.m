function s = getConstraint(x,y)
global object radius

if object == "Ball"
    if x > -radius
        s = sqrt(radius^2 - x.^2);
    elseif x <= -radius
        s = 0;
    end
elseif object == "Cubic"
    if x > -radius
        s = radius * (y >= -0.18 && y <= -0.12);
    elseif x <= -radius
        s = 0;
    end
end
end