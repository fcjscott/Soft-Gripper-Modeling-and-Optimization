syms xkm1 ykm1 xk yk xkp1 ykp1 tau phi_k0 dl EI kappa F k



e0 = [xk-xkm1;yk-ykm1];
e1 = [xkp1-xk;ykp1-yk];
e0Dote1 = (xk-xkm1)*(xkp1-xk)+(yk-ykm1)*(ykp1-yk);
e0Crosse1 = e1(1)*e0(2)-e1(2)*e0(1);
phi_k = atan(e0Crosse1/e0Dote1);
curvature = 2*tan(phi_k/2);
curvature_0 = kappa;       

F1 = diff((curvature - curvature_0)^2,xkm1)
F2 = diff((curvature - curvature_0)^2,ykm1)
F3 = diff((curvature - curvature_0)^2,xk)
F4 = diff((curvature - curvature_0)^2,yk)
F5 = diff((curvature - curvature_0)^2,xkp1)
F6 = diff((curvature - curvature_0)^2,ykp1)
%%
J11 = diff(F1,xkm1)
J12 = diff(F1,ykm1)
J13 = diff(F1,xk)
J14 = diff(F1,yk)
J15 = diff(F1,xkp1)
J16 = diff(F1,ykp1)
J22 = diff(F2,ykm1)
J23 = diff(F2,xk)
J24 = diff(F2,yk)
J25 = diff(F2,xkp1)
J26 = diff(F2,ykp1)
J33 = diff(F3,xk)
J34 = diff(F3,yk)
J35 = diff(F3,xkp1)
J36 = diff(F3,ykp1)
J44 = diff(F4,yk)
J45 = diff(F4,xkp1)
J46 = diff(F4,ykp1)
J55 = diff(F5,xkp1)
J56 = diff(F5,ykp1)
J66 = diff(F6,ykp1)


