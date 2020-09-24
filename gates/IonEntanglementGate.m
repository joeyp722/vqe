function RhoOut = IonEntanglementGate(RhoIn,param)
%Executes sequence of entanglement gates with or without decoherence
%based on the error flag based on the ion implementation, here a molmer
%sorensen gate is implemented in the X X Pauli bases in chain form.
%input: RhoIn, param, angles
%these are:
%input density matrix, parameters related to implementation type
%output: RhoOut
%this is: output density matrix

v_plus=[1; 1]/sqrt(2);
v_min=[1; -1]/sqrt(2);

GeometricPhaseGate=kron(v_plus,v_min)*kron(v_plus,v_min)'+kron(v_min,v_plus)*kron(v_min,v_plus)'-1i*kron(v_plus,v_plus)*kron(v_plus,v_plus)'-1i*kron(v_min,v_min)*kron(v_min,v_min)';

Rho=RhoIn;

nqbits=log(length(Rho))/log(2);

if nqbits==2
    Rho=GeometricPhaseGate*Rho*GeometricPhaseGate';
else
    for i=1:nqbits-1
        M = kron(eye(2^(i-1)),kron(GeometricPhaseGate,eye(2^(nqbits-i-1))));
        Rho = M*Rho*M';
        if param.error; Rho=LindbladSelect(Rho,param,'entanglement'); end
    end
end

RhoOut=Rho;
