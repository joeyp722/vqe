function output = vqe_angles(param)
%vqe subroutine that solves the potential and fidelity for the ground state
%of a molecule
%param is input struct
%inputs are:
%step, depth, name, distance, hamiltonian, initial angles
%these are the:
%number of iterations, circuit depth, name of the molecule, distance between atoms, hamiltonian in the Pauli
%basis
%output is the output struct
%outputs are:
%Potential, Energy, Rho, Fidelity, kvec, potential_curve, pseudo_curve,
%Pseudo_potential
%these are the:
%Energy potential, Energy potential without repulsion and ionization
%energies, density matrix of the ground state, fldelity between the this subroutine and matlabs eigensolver, 
%axis of iterations, potential with respect to kvec, 
%pseudo potential with respect to kvec, potential calculated by the maltab eigensolver

%Setting parameters with param input struct
steps = param.steps; %Number of iterations of algorithm
depth = param.depth; %Circuit depth
name = param.name; %Name of molecule

distance = param.distance; %Distance between atoms in molecule

H=param.hamiltonian; %Hamiltonian that is solved

nqbits=log(length(H))/log(2); %Number of required qubits based on Hamiltonian
t=param.angles; %Initial angles

%Initial density matrix
RhoInit=zeros(2^nqbits,2^nqbits);
RhoInit(end,end)=1;

%SPSA parameters
alpha=param.alpha;
gamma=param.gamma;
c=param.c;

%Calculation of additional SPSA parameter
average_start_grad=get_average_start_grad(param);
a=2*pi*c/(5*average_start_grad);

%Start of the algorithm
for k=1:1:steps

deltak=randi([0 1],depth+1,3*nqbits); %Random direction increment
ck=c/(k^gamma); %Updating SPSA parameter

angles=rem(real(t+ck*deltak),2*pi); %Incrementing angles in the positive direction

%Quantum circuit
Rho=RhoInit; %Reseting density matrix
Rho=MultiRotationGateSelect(Rho,param,angles(1,:)); %Single qubit rotational sequence

for i=2:depth+1
    Rho=EntanglementGateSelect(Rho,param); %Entangler gate sequence
    Rho=MultiRotationGateSelect(Rho,param,angles(i,:)); %Single qubit rotational sequence
end

%Measurement of cost function
Energyplus=trace(Rho*H); 

angles=rem(real(t-ck*deltak),2*pi); %Incrementing angles in the negative direction

%Quantum circuit
Rho=RhoInit; %Reseting density matrix
Rho=MultiRotationGateSelect(Rho,param,angles(1,:)); %Single qubit rotational sequence

for i=2:depth+1
    Rho=EntanglementGateSelect(Rho,param); %Entangler gate sequence
    Rho=MultiRotationGateSelect(Rho,param,angles(i,:)); %Single qubit rotational sequence
end

Energymin=trace(Rho*H); %Measurement of cost function

gk=(Energyplus-Energymin)/(2*ck)*deltak; %Updating SPSA parameter based on cost functions
ak=a/(k^alpha); %Updating SPSA parameter
Energy=0.5*(Energyplus+Energymin); %Solution based on cost function
energy(k,1)=Energy; %Solution as function of iteration
kvec(k,1)=k; %Iteration axis

t=t-ak*gk; %Update angles
end 
Energy=0.5*(Energyplus+Energymin); %Solution based on cost function
Repulsion=get_Repulsion(name,distance);
Ionization_Energies=get_Ionization_Energies(name);
Potential=(real(Energy)+Repulsion+Ionization_Energies); %Calculating potential
potential_curve=(energy+(Repulsion+Ionization_Energies)*ones(steps,1)); %Calculating potential as function of iteration

%Calculate eigenvalues and eigenvectors
[eigen_vec,eigen_energy]=eig(H);
[eigen_energy,index] = sort(diag(eigen_energy));
eigen_vec = eigen_vec(:, index);
Rho_nonvqe=eigen_vec(:,1)*eigen_vec(:,1)';

%Calculating pseudo potential
pseudo_Energy=eigen_energy(1);
pseudo_Potential=(real(pseudo_Energy)+Repulsion+Ionization_Energies); %Calculating pseudo potential

%Constructing output struct
output.kvec=kvec;
output.potential_curve=potential_curve;
output.Potential=Potential;
output.pseudo_Potential=pseudo_Potential;

output.Energy=Energy;
output.Rho=Rho;
output.Fidelity=get_fidelity(Rho,Rho_nonvqe); %Calculating fidelity