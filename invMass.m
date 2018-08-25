function im = invMass(x,y)
im = eye(2,2) - (nConstraint(x,y))*(nConstraint(x,y))';
end
