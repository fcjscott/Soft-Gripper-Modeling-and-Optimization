function n = nConstraint(x,y)
global object radius

if object == "Ball"
    n = [x;y + 0.15]/radius;
elseif object == "Cubic"
    n = [-1;0];
end

end
