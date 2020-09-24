function RhoOut = MultiTransmonRotationGate(RhoIn,param,angles)
%Executes sequence of ZXZ rotational gates with or without decoherence
%based on the error flag based on the transmon implementation.
%input: RhoIn, param, angles
%these are:
%input density matrix, parameters related to implementation type, rotational angles of
%the gates (this is a matrix)
%output: RhoOut
%this is: output density matrix

Rho=RhoIn;

nqbits=length(angles)/3;
angles = reshape(angles,[nqbits,3]);

rabi=param.rabi;

I=[1 0;0 1];
X=[0 1;1 0];
Y=[0 -1i;1i 0];
Z=[1 0;0 -1];


U=1;

for i=1:nqbits
    U_z1=cos(angles(i,1)/2)*I-1i*sin(angles(i,1)/2)*Z;
    U=kron(U,U_z1);
end

Rho=U*Rho*U';
param.rotational_time=max(abs(real(angles(:,1))))/rabi;
if param.error; Rho=LindbladSelect(Rho,param,'rotational'); end


U=1;

for i=1:nqbits
    U_x=cos(angles(i,2)/2)*I-1i*sin(angles(i,2)/2)*X;
    U=kron(U,U_x);
end

Rho=U*Rho*U';
param.rotational_time=max(abs(real(angles(:,2))))/rabi;
if param.error; Rho=LindbladSelect(Rho,param,'rotational'); end


U=1;

for i=1:nqbits
    U_z2=cos(angles(i,3)/2)*I-1i*sin(angles(i,3)/2)*Z;
    U=kron(U,U_z2);
end

Rho=U*Rho*U';
param.rotational_time=max(abs(real(angles(:,3))))/rabi;
if param.error; Rho=LindbladSelect(Rho,param,'rotational'); end

RhoOut=Rho;