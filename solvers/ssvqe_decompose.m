function output = ssvqe_decompose(param)
%ssvqe (subspace vqe) subroutine that solves the potential and fidelity for all states
%using the ssvqe weighted method where the entanglers are also optimized
%of a molecule
%param is input struct
%inputs are:
%step, depth, name, distance, hamiltonian
%these are the:
%number of iterations, circuit depth, name of the molecule, distance between atoms, hamiltonian in the Pauli
%basis
%output is the output struct
%outputs are:
%Potential, Energy, Rho, Fidelity, Pseudo_potential
%these are the:
%Energy potential, Energy potential without repulsion and ionization
%energies, density matrix of the ground state, fldelity between the this subroutine and matlabs eigensolver,  
%potential calculated by the maltab eigensolver

%Setting parameters with param input struct
steps = param.steps; %Number of iterations of algorithm
depth = param.depth; %Circuit depth
name = param.name; %Name of molecule

distance = param.distance; %Distance between atoms in molecule

H=param.hamiltonian; %Hamiltonian that is solved

nstates=length(H); %Number of states  based on Hamiltonian
nqbits=log(length(H))/log(2); %Number of required qubits based on Hamiltonian
t=ones(2*depth+1,3*nqbits); %Initial angles

weight_multiplier=param.weight_multiplier; %Mulitplier used within cost function

%SPSA parameters
alpha=param.alpha;
gamma=param.gamma;
c=param.c;

%Calculation of additional SPSA parameter
average_start_grad=get_average_start_grad(param);
a=2*pi*c/(5*average_start_grad);

%Start of the algorithm
for k=1:1:steps
    
deltak=randi([0 1],2*depth+1,3*nqbits); %Random direction increment
ck=c/(k^gamma); %Updating SPSA paramete

angles=t+ck*deltak; %Incrementing angles in the positive direction

Energyplus=zeros(1,nstates); %Initializing array for storing calculated energies for cost function

%Evaluation of every state
for j=1:nstates
    RhoInit=zeros(2^nqbits,2^nqbits); 
    RhoInit(j,j)=1; %Initial density matrix for the j-th state
    
    %Quantum circuit
    Rho=RhoInit; %Reseting density matrix
    Rho=MultiRotationGateSelect(Rho,param,angles(1,:)); %Single qubit rotational sequence

    for i=1:depth
        Rho=EntanglementRotationGate(Rho,param,angles(2*i,:)); %Entangler gate sequence
        Rho=MultiRotationGateSelect(Rho,param,angles(2*i+1,:)); %Single qubit rotational sequence
    end
    
    Energyplus(j)=trace(Rho*H); %Storing calculated energies for cost function
end

angles=t-ck*deltak; %Incrementing angles in the negative direction

Energymin=zeros(1,nstates); %Initializing array for storing calculated energies for cost function

%Evaluation of every state
for j=1:nstates
    RhoInit=zeros(2^nqbits,2^nqbits);
    RhoInit(j,j)=1; %Initial density matrix for the j-th state
    
    %Quantum circuit
    Rho=RhoInit; %Reseting density matrix
    Rho=MultiRotationGateSelect(Rho,param,angles(1,:)); %Single qubit rotational sequence

    for i=1:depth
        Rho=EntanglementRotationGate(Rho,param,angles(2*i,:)); %Entangler gate sequence
        Rho=MultiRotationGateSelect(Rho,param,angles(2*i+1,:)); %Single qubit rotational sequence
    end
    
    %Storing the density matrix for the state being evaluated
    if k==steps
        data.Rho{j}=Rho;
    end
    
    Energymin(j)=trace(Rho*H); %Storing calculated energies for cost function
end

%Evaluation of cost function
Energymin_weight=0; %Setting variable for cost function where angles increment in the negative direction
Energyplus_weight=0; %Setting variable for cost function where angles increment in the positive direction

Energy=zeros(1,nstates); %Initializing storage for solution of all the states
for j=1:nstates
    w=j;
    Energy(j)=0.5*(Energyplus(j)+Energymin(j)); %Solution of the j-th state
    Energymin_weight=Energymin_weight+weight_multiplier*w*Energymin(j); %Cost function where angles increment in the negative direction
    Energyplus_weight=Energyplus_weight+weight_multiplier*w*Energyplus(j); %Cost function where angles increment in the positive direction  
end

gk=(Energyplus_weight-Energymin_weight)/(2*ck)*deltak; %Updating SPSA parameter based on cost functions
ak=a/(k^alpha); %Updating SPSA parameter

t=t-ak*gk; %Update angles

end 

Repulsion=get_Repulsion(name,distance);
Ionization_Energies=get_Ionization_Energies(name);
Potential=(real(Energy)+(Repulsion+Ionization_Energies)*ones(1,length(Energy))); %Calculating potentials

[Potential,index]=sort(Potential); %Sorting Potentials based, where the ordering is saved in index

for j=1:nstates
    output.Rho{j}=data.Rho{index(j)}; %Constructing output.Rho struct where the density matrices are sorted based on index
end

%Calculate eigenvalues and eigenvectors
[eigen_vec,eigen_energy]=eig(H);
[eigen_energy,index] = sort(diag(eigen_energy));
eigen_vec = eigen_vec(:, index);
    
%Calculating pseudo potentials
for j=1:nstates
    output.nonvqe_Rho{j}=eigen_vec(:,j)*eigen_vec(:,j)';
    pseudo_Energy(j)=eigen_energy(j);
    pseudo_Potential(j)=(real(pseudo_Energy(j))+(Repulsion+Ionization_Energies)*ones(1,length(pseudo_Energy(j))));
end

%Constructing output struct
for j=1:nstates
    output.Potential{j}=Potential(j);
    output.pseudo_Potential{j}=pseudo_Potential(j);
    output.Fidelity{j}=get_fidelity(output.nonvqe_Rho{j},output.Rho{j});
    get_fidelity(output.nonvqe_Rho{j},output.Rho{j}) %Calculating fidelity
end