function Fy = getLift
global ForceAll object miu nv q

fNormal = zeros(nv-2,1);
Fup = zeros(nv-2,1);
for c = 3:nv
    xPos = q(2*c-1);
    yPos = q(2*c);
    if strcmp(object,"Cubic")
        Fup(c-2) = - (ForceAll(2*c-1)*miu + ForceAll(2*c));
    elseif strcmp(object,"Ball")
        fReaction = ForceAll(2*c-1:2*c);
        fNormal(c-2) = dot(fReaction, nConstraint(xPos,yPos));
        Fup(c-2) = - ForceAll(2*c) + fNormal(c-2)*miu*cosd...
        (atand(abs(yPos+0.15)/abs(xPos)));
    end
end
Fy = sum(Fup);        

end