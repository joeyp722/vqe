function RhoOut = MultiMajoranaRotationGate(RhoIn,param,angles)
%Executes sequence of XYX rotational gates with or without decoherence
%based on the error flag based of the majorana implementation.
%input: RhoIn, param, angles
%these are:
%input density matrix, parameters related to implementation type, rotational angles of
%the gates (this is a matrix)
%output: RhoOut
%this is: output density matrix

Rho=RhoIn;

nqbits=length(angles)/3;
angles = reshape(angles,[nqbits,3]);

I=[1 0;0 1];
X=[0 1;1 0];
Y=[0 -1i;1i 0];
Z=[1 0;0 -1];


U=1;

for i=1:nqbits
    U_x1=cos(angles(i,1)/2)*I-1i*sin(angles(i,1)/2)*X;
    U=kron(U,U_x1);
end

Rho=U*Rho*U';
if param.error; Rho=LindbladSelect(Rho,param,'rotational'); end


U=1;

for i=1:nqbits
    U_y=cos(angles(i,2)/2)*I-1i*sin(angles(i,2)/2)*Y;
    U=kron(U,U_y);
end

Rho=U*Rho*U';
if param.error; Rho=LindbladSelect(Rho,param,'rotational'); end


U=1;

for i=1:nqbits
    U_x2=cos(angles(i,3)/2)*I-1i*sin(angles(i,3)/2)*X;
    U=kron(U,U_x2);
end

Rho=U*Rho*U';
if param.error; Rho=LindbladSelect(Rho,param,'rotational'); end

RhoOut=Rho;