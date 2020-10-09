function [Kg,Mg]=element_m(E,rho,L,A,I,R,Ax,Ix)
%mass and stiffness matrices of linear beam element
%transformation is given by R = [cos,sin]

%sensitivities
if nargin>6, A = Ax; I = Ix; end

%transformation matrix
T = transformation(R);

%Euler-Bernoulli beam element
AEL = A*E/L; EIL1 = 2*E*I/L; EIL2 = 6*E*I/L^2; EIL3 = 12*E*I/L^3; EIL4 = 4*E*I/L;

Ke = [
    AEL,     0,     0, -AEL,     0,     0
    0,    EIL3,  EIL2,    0, -EIL3,  EIL2
    0,    EIL2,  EIL4,    0, -EIL2,  EIL1
    -AEL,    0,     0,  AEL,     0,     0
    0,   -EIL3, -EIL2,    0,  EIL3, -EIL2
    0,    EIL2,  EIL1,    0, -EIL2,  EIL4 ];


%transformation of coordinates
Kg = T.'*Ke*T;

Me = rho*A/420*[
    140*L,     0,      0,  70*L,       0,       0
    0,     156*L, 22*L^2,     0,    54*L, -13*L^2
    0,    22*L^2,  4*L^3,     0,  13*L^2,  -3*L^3
    70*L,      0,      0, 140*L,       0,       0
    0,      54*L, 13*L^2,     0,   156*L, -22*L^2
    0,   -13*L^2, -3*L^3,     0, -22*L^2,   4*L^3 ];

%transformation of coordinates
Mg = T.'*Me*T;
    
end
