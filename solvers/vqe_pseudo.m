function output = vqe_pseudo(param)
%pseudo vqe subroutine that uses the matlab eigensolver
%param is input struct
%inputs are:
%name, distance, hamiltonian
%these are the:
%name of the molecule, distance between atoms, hamiltonian in the Pauli
%basis
%output is the output struct
%outputs are:
%Potential, Energy, Rho, Fidelity
%these are the:
%Energy potential, Energy potential without repulsion and ionization
%energies, density matrix of ground state, fldelity of 100%

%Setting parameters with param input struct
name = param.name;
distance = param.distance;

H=param.hamiltonian;

%Calculate eigenvalues and eigenvectors
[eigen_vec,eigen_energy]=eig(H);
[eigen_energy,index] = sort(diag(eigen_energy));
eigen_vec = eigen_vec(:, index);
Rho=eigen_vec(:,1)*eigen_vec(:,1)'; %Constructing density matrix ground state

%Calcualting potential
Energy=eigen_energy(1); %Setting ground state energy
Repulsion=get_Repulsion(name,distance);
Ionization_Energies=get_Ionization_Energies(name);
Potential=(real(Energy)+Repulsion+Ionization_Energies); %Calculating potential

%Constructing output struct
output.Potential=Potential;
output.Energy=Energy;
output.Rho=Rho;
output.Fidelity=1;