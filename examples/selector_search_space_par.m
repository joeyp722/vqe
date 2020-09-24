%example script that plots the potential of molecular hydrogen near its
%bond length and the fidelity with respect to the matlab eigensolver
%solution as a result of the initial angles of the vqe, the angles are chosen at random. This
%example uses parallel processing.

clear all;

%select implemtation type
select='General';
param = LoadParam(select);
param.select=select;

%set iterations for repition of the vqe subroutine
iterations=5;

%set vqe subroute parameters
param.steps=2000; %iterations until solution is found
param.depth=3; %circuit depth

%SPSA parameters
param.alpha=0.602;
param.gamma=0.101;
param.c=0.1;

%parameters for construct hamiltonian
param.name='H2'; %molecule name
param.transformation_type='jw'; %hamiltonian construction method
param.electron_number=2; %number of valence electrons
param.spin_number=0; %total valence spin number
param.swap=true; %swap ground and excited states in the hamiltonian

%error flag
param.error=false; %error flag

%ditance parameters
param.distance=0.75; %set distance
param.distance_interval=0.05;
param.min_distance=0.10;
param.max_distance=2.30;

distance_interval=param.distance_interval;
min_distance=param.min_distance;
max_distance=param.max_distance;

distance_index=(param.distance+distance_interval-min_distance)/distance_interval;

name=param.name;

%load matrices
load(strcat('molecules/',name,'/one_electron_matrix',num2str(distance_index),'.mat'));
one_electron_matrix=matrix;
load(strcat('molecules/',name,'/two_electron_matrix',num2str(distance_index),'.mat'));
two_electron_matrix=matrix;

%construct hamiltonian
 param.hamiltonian=get_hamiltonian(one_electron_matrix,two_electron_matrix,param);
 if param.swap; param.hamiltonian=swap(param.hamiltonian); end

%set counter variables
count_array=[length(relax_exp),length(dephase_exp),iterations];
count_total=prod(count_array);

%iterate for different decay rates
parfor_progress(iterations);
parfor i=1:iterations
    
    %Constructing param for each parallel process
    param_par(i)=param;
    %vqe subroutine
    
    %Set random initial angles
    nqbits=log(length(param_par(i).hamiltonian))/log(2);
    param_par(i).angles=2*pi*rand(param_par(i).depth+1,3*nqbits)
    
    output_par(i)=vqe_angles(param_par(i));
    parfor_progress;
end
parfor_progress(0);

%Unpack parallel processing results
for i=1:iterations
     
    output{i}=output_par(i);
    
end

axis_iter=[];
for i=1:iterations
	axis_iter(i)=i;
end

%calculating the mean and std of the potential based on multiple iterations
dataset=[];
for i=1:iterations
	dataset(i)=output{i}.Potential;
end
output{1}.potential.mean=mean(dataset);
output{1}.potential.std=std(dataset);
potential_dataset=dataset;

%calculating the mean and std of the fidelity based on multiple iterations
dataset=[];
for i=1:iterations
    dataset(i)=output{i}.Fidelity;
end
output{1}.fidelity.mean=mean(dataset);
output{1}.fidelity.std=std(dataset);
fidelity_dataset=dataset;

Potential_string=strcat('Potential mean: ',num2str(mean(potential_dataset)),', Potential std: ',num2str(std(potential_dataset)));
Fidelity_string=strcat('Fidelity mean: ',num2str(mean(fidelity_dataset)),', Fidelity std: ',num2str(std(fidelity_dataset)));
disp(Potential_string);
disp(Fidelity_string)

%plotting results
figure(1);
scatter(axis_iter,potential_dataset);
title('Potential');
xlabel('Iteration');
ylabel('Potential');
legend(Potential_string);
savefig('Potential_scatter.fig');


figure(2);
scatter(axis_iter,fidelity_dataset);
title('Fidelity');
xlabel('Iteration');
ylabel('Fidelity');
legend(Fidelity_string);
savefig('Fidelity_scatter.fig');