function output = get_average_start_grad(param)
%Gives the average energy solution for the vqe algorithm in the case of 25
%random directions.
%input: param struct: depth, name, hamiltonian
%these are: circuit depth, name of the molecule, Hamiltonian of molecule
%output: output
%this is: stated above

depth = param.depth;
name = param.name;

H=param.hamiltonian;

nqbits=log(length(H))/log(2);
t=ones(depth+1,3*nqbits);

RhoInit=zeros(2^nqbits,2^nqbits);
RhoInit(end,end)=1;

c=0.1;

for k=1:1:25
deltak=randi([0 1],depth+1,3*nqbits);
angles=rem(real(t+c*deltak),2*pi);

Rho=RhoInit;
Rho=MultiRotationGateSelect(Rho,param,angles(1,:));

for i=2:depth+1
    Rho=EntanglementGateSelect(Rho,param);
    Rho=MultiRotationGateSelect(Rho,param,angles(i,:));
end

Energyplus=trace(Rho*H);


angles=rem(real(t-c*deltak),2*pi);

Rho=RhoInit;
Rho=MultiRotationGateSelect(Rho,param,angles(1,:));

for i=2:depth+1
    Rho=EntanglementGateSelect(Rho,param);  
    Rho=MultiRotationGateSelect(Rho,param,angles(i,:));
end

Energymin=trace(Rho*H);

Energy=(Energyplus-Energymin);
energy(k,1)=Energy;

end 

output=mean(abs(real(energy)));
