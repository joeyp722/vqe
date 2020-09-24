function RhoOut = MajoranaEntanglementGate(RhoIn,param)
%Executes sequence of entanglement gates with or without decoherence
%based on the error flag based on the majorana implementation, here a cnot 
%gate is implemented in chain form.
%input: RhoIn, param, angles
%these are:
%input density matrix, parameters related to implementation type
%output: RhoOut
%this is: output density matrix

EntanglementGate=[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0];

Rho=RhoIn;

nqbits=log(length(Rho))/log(2);

if nqbits==2
    Rho=EntanglementGate*Rho*EntanglementGate';
else
    for i=1:nqbits-1
        M = kron(eye(2^(i-1)),kron(EntanglementGate,eye(2^(nqbits-i-1))));
        Rho = M*Rho*M';
        if param.error; Rho=LindbladSelect(Rho,param,'entanglement'); end
    end
end

RhoOut=Rho;
