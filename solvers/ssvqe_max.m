function output = ssvqe_max(param)
%ssvqe (subspace vqe) subroutine that solves the potential and fidelity for all states
%using the ssvqe max method
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
%potential calculated by the maltab eigensolver%

%Setting parameters with param input struct
steps = param.steps; %Number of iterations of algorithm
depth = param.depth; %Circuit depth
name = param.name; %Name of molecule

distance = param.distance; %Distance between atoms in molecule
 
H=param.hamiltonian; %Hamiltonian that is solved

filter_weight=param.filter_weight; %Mulitplier used within cost function

nstates=length(H); %Number of states  based on Hamiltonian
nqbits=log(length(H))/log(2); %Number of required qubits based on Hamiltonian
t=ones(depth+1,3*nqbits); %Initial angles

Energy=zeros(1,nstates); %Initializing storage for solution of all the states

%SPSA parameters
alpha=param.alpha;
gamma=param.gamma;
c=param.c;

%Calculation of additional SPSA parameter
average_start_grad=get_average_start_grad(param);
a=2*pi*c/(5*average_start_grad);

%Evaluation of every state
for j=1:nstates
    
    %First part algorithm
    
    %Constructing initial density matrix
    RhoInit=zeros(2^nqbits,2^nqbits);
    for i=1:j
        RhoInit(i,i)=1/j;
    end
    
    %Initialization of filter for the cost function
    Filter=zeros(2^nqbits,2^nqbits);
    for i=(j+1):nstates
        Filter(i,i)=1;
    end

    %Start of the first part of the algorithm
    for k=1:1:steps
    deltak=randi([0 1],depth+1,3*nqbits); %Random direction increment
    ck=c/(k^gamma); %Updating SPSA parameter

    anglesu=rem(real(t+ck*deltak),2*pi); %Incrementing angles in the positive direction

    %Quantum circuit
    Rho=RhoInit; %Reseting density matrix
    Rho=MultiRotationGateSelect(Rho,param,anglesu(1,:)); %Single qubit rotational sequence

    for i=2:depth+1
        Rho=EntanglementGateSelect(Rho,param); %Entangler gate sequence
        Rho=MultiRotationGateSelect(Rho,param,anglesu(i,:)); %Single qubit rotational sequence
    end
    
    %Measurement of cost function
    Energyplus=trace(Rho*H);

    anglesu=rem(real(t-ck*deltak),2*pi); %Incrementing angles in the negative direction

    %Quantum circuit
    Rho=RhoInit; %Reseting density matrix
    Rho=MultiRotationGateSelect(Rho,param,anglesu(1,:)); %Single qubit rotational sequence

    for i=2:depth+1
        Rho=EntanglementGateSelect(Rho,param); %Entangler gate sequenc
        Rho=MultiRotationGateSelect(Rho,param,anglesu(i,:)); %Single qubit rotational sequence
    end
    
    %Measurement of cost function
    Energymin=trace(Rho*H);

    gk=(Energyplus-Energymin)/(2*ck)*deltak; %Updating SPSA parameter based on cost functions
    ak=a/(k^alpha); %Updating SPSA parameter
    
    t=t-ak*gk; %Update angles
    end

    %Second part algorithm
    
    t=ones(depth+1,3*nqbits); %Initial angles
    
    %Constructing initial density matrix
    RhoInit=zeros(2^nqbits,2^nqbits);
    RhoInit(end,end)=1;
    
    if j ~=1
        
        %Initialization of filter for the cost function
        Filter=zeros(2^nqbits,2^nqbits);
        for i=(j+1):nstates
            Filter(i,i)=1;
        end
        
        %Start of the second part of the algorithm
        for k=1:1:steps
        deltak=randi([0 1],depth+1,3*nqbits); %Random direction increment
        ck=c/(k^gamma); %Updating SPSA parameter

        
        angless=rem(real(t+ck*deltak),2*pi); %Incrementing angles in the positive direction

        %Quantum circuit
        Rho=RhoInit; %Reseting density matrix
        Rho=MultiRotationGateSelect(Rho,param,angless(1,:)); %Single qubit rotational sequence

        for i=2:depth+1
            Rho=EntanglementGateSelect(Rho,param); %Entangler gate sequence
            Rho=MultiRotationGateSelect(Rho,param,angless(i,:)); %Single qubit rotational sequence
        end
        
        %Measurement of the filter part of the cost function
        CostSubspaceplus=trace(Rho*Filter);
        Rho=MultiRotationGateSelect(Rho,param,anglesu(1,:)); %Single qubit rotational sequence

        for i=2:depth+1
            Rho=EntanglementGateSelect(Rho,param); %Entangler gate sequence
            Rho=MultiRotationGateSelect(Rho,param,anglesu(i,:)); %Single qubit rotational sequence
        end

        Energyplus=trace(Rho*H); %Measurement of the energy part of the cost function

        angless=rem(real(t-ck*deltak),2*pi); %Incrementing angles in the negative direction

        %Quantum circuit
        Rho=RhoInit; %Reseting density matrix
        Rho=MultiRotationGateSelect(Rho,param,angless(1,:)); %Single qubit rotational sequence

        for i=2:depth+1
            Rho=EntanglementGateSelect(Rho,param); %Entangler gate sequence
            Rho=MultiRotationGateSelect(Rho,param,angless(i,:)); %Single qubit rotational sequence
        end

        %Measurement of the filter part of the cost function
        CostSubspacemin=trace(Rho*Filter);
        Rho=MultiRotationGateSelect(Rho,param,anglesu(1,:)); %Single qubit rotational sequence

        for i=2:depth+1
            Rho=EntanglementGateSelect(Rho,param); %Entangler gate sequence
            Rho=MultiRotationGateSelect(Rho,param,anglesu(i,:)); %Single qubit rotational sequence
        end

        Energymin=trace(Rho*H); %Measurement of the energy part of the cost function

        %Storing the density matrix for the state being evaluated
        if k==steps
            data.Rho{j}=Rho;
        end
        
        %Evaluation of cost function
        CostEnergy=Energyplus-Energymin;
        CostSubspace=CostSubspaceplus-CostSubspacemin;

        gk=(CostSubspace*filter_weight-CostEnergy)/(2*ck)*deltak; %Updating SPSA parameter based on cost functions
        ak=a/(k^alpha); %Updating SPSA parameter

        t=t-ak*gk; %Update angles
        end
    end
    data.Rho{j}=Rho; %Storing the density matrix for the state being evaluated
    Energy(j)=0.5*(Energyplus+Energymin); %Solution of the j-th state
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
    get_fidelity(output.nonvqe_Rho{j},output.Rho{j})
end