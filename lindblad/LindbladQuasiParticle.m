function RhoOut = LindbladQuasiParticle(RhoIn,p,p_qp,p_cor_even,p_cor_odd,p_pair)
%Kraus transformation for the quasi particle poisoning channel for the majorana
%of the qubit.
%input: RhoIn, ... , probabilities
%these are: input density matrix, ... , probablities related to the
%poisoning
%ouput: RhoOut
%this is: output density matrix

Rho=RhoIn;

nqbits=log(length(Rho))/log(2);
p_ne=1-p;

E0=sqrt(p_ne+p_qp+p_cor_even)*[1 0;0 1];
E1=sqrt(p_cor_odd)*[0 1;1 0];
E2=sqrt(p_pair)*[1 0;0 -1];

for i=1:nqbits
    M0=kron(eye(2^(i-1)),kron(E0,eye(2^(nqbits-i))));
    M1=kron(eye(2^(i-1)),kron(E1,eye(2^(nqbits-i))));
    M2=kron(eye(2^(i-1)),kron(E2,eye(2^(nqbits-i))));
    Rho=M0*Rho*M0'+M1*Rho*M1'+M2*Rho*M2';
end

RhoOut=Rho;