function output = get_fidelity(rho1,rho2)
%Gives fidelity between rho1 and rho2
%input: rho1, rho2
%these are: density matrix 1, density matrix 2
%output: output
%this is: the Uhlmann fidelity

offset=10^-10; %To prevent: Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.
fidelity=real(trace(((rho1+offset)^0.5*rho2*(rho1+offset)^0.5)^0.5))^2;

output=fidelity;