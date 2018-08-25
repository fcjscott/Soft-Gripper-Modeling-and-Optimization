function [P,Ps] = getPressure(t,Nsteps,tempP)
global dt Pd

Tfinal = t*dt;
damp = 0.707;
ts = Nsteps/2*dt;
omegaN = 4/(ts*damp);
nume = omegaN^2;
deno = [1 2*damp*omegaN omegaN^2];
sys = tf(nume,deno);
Ps = lsim(sys,[tempP;Pd],0:dt:Tfinal);
P = Ps(end);

end
