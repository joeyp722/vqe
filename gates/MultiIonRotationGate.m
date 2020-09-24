function RhoOut = MultiIonRotationGate(RhoIn,param,angles)
%Executes sequence of rotational gates with or without decoherence
%based on the error flag based on the ion implementation
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
detuning=param.detuning;
rabi_eff=sqrt(rabi^2+detuning^2);

I=[1 0;0 1];
X=[0 1;1 0];
Y=[0 -1i;1i 0];
Z=[1 0;0 -1];


U=1;

for i=1:nqbits
    U_theta1=cos(angles(i,1)/2)*I-1i*sin(angles(i,1)/2)*(rabi*X-detuning*Z)/rabi_eff;
    U=kron(U,U_theta1);
end

Rho=U*Rho*U';
param.rotational_time=max(abs(real(angles(:,1))))/rabi;
if param.error; Rho=LindbladSelect(Rho,param,'rotational'); end


U=1;

for i=1:nqbits
    U_r=cos(angles(i,2)/2)*I-1i*sin(angles(i,2)/2)*Z;
    U=kron(U,U_r);
end

Rho=U*Rho*U';
param.rotational_time=max(abs(real(angles(:,2))))/rabi;
if param.error; Rho=LindbladSelect(Rho,param,'rotational'); end


U=1;

for i=1:nqbits
    U_theta2=cos(angles(i,3)/2)*I-1i*sin(angles(i,3)/2)*(rabi*X-detuning*Z)/rabi_eff;
    U=kron(U,U_theta2);
end

Rho=U*Rho*U';
param.rotational_time=max(abs(real(angles(:,3))))/rabi;
if param.error; Rho=LindbladSelect(Rho,param,'rotational'); end

RhoOut=Rho;