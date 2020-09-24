%example script that plots the potential of molecular hydrogen near its
%bond length and the fidelity with respect to the matlab eigensolver
%solution as a function of different relaxation and dephasing rates.

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

%error flag and decay rates
param.error=true; %error flag
relax_exp=-6:1:0; %relaxation rate
dephase_exp=-4:1:6; %dephasing rate

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

%iterate for different decay rates
for relax_index=1:length(relax_exp)
    for dephase_index=1:length(dephase_exp)
        for i=1:iterations
            param.relax=10^relax_exp(relax_index);
            param.dephase=10^dephase_exp(dephase_index);
            %vqe subroutine
            output{i,relax_index,dephase_index}=vqe(param);
        end
        disp([relax_exp(relax_index),dephase_exp(dephase_index),output{i,relax_index,dephase_index}.Fidelity]);
    end
end

%calculating the mean and std of the potential based on multiple iterations
for relax_index=1:length(relax_exp)
    for dephase_index=1:length(dephase_exp)
        dataset=[];
        for i=1:iterations
            dataset(i)=output{i,relax_index,dephase_index}.Potential;
        end
        output{1,relax_index,dephase_index}.potential.mean=mean(dataset);
        output{1,relax_index,dephase_index}.potential.std=std(dataset);
    end
end

%calculating the mean and std of the fidelity based on multiple iterations
for relax_index=1:length(relax_exp)
    for dephase_index=1:length(dephase_exp)
        dataset=[];
        for i=1:iterations
            dataset(i)=output{i,relax_index,dephase_index}.Fidelity;
        end
        output{1,relax_index,dephase_index}.fidelity.mean=mean(dataset);
        output{1,relax_index,dephase_index}.fidelity.std=std(dataset);
    end
end

%defining plotting axes
potential_mean=[];
potential_std=[];
fidelity_mean=[];
fidelity_std=[];

%Putting the data in the plotting axes
for relax_index=1:length(relax_exp)
    for dephase_index=1:length(dephase_exp)
        potential_mean(relax_index,dephase_index)=output{1,relax_index,dephase_index}.potential.mean;
        potential_std(relax_index,dephase_index)=output{1,relax_index,dephase_index}.potential.std;
        fidelity_mean(relax_index,dephase_index)=output{1,relax_index,dephase_index}.fidelity.mean;
        fidelity_std(relax_index,dephase_index)=output{1,relax_index,dephase_index}.fidelity.std;
    end
end


%plotting results
figure(1);
surf(dephase_exp,relax_exp,potential_mean);
title('Potential mean');
xlabel('log_{10}(dephase)');
ylabel('log_{10}(relax)');
zlabel('Potential [Hartrees]');

figure(2);
surf(dephase_exp,relax_exp,potential_std);
title('Potential std');
xlabel('log_{10}(dephase)');
ylabel('log_{10}(relax)');
zlabel('Potential [Hartrees]');

figure(3);
surf(dephase_exp,relax_exp,fidelity_mean);
title('Fidelity mean');
xlabel('log_{10}(dephase)');
ylabel('log_{10}(relax)');
zlabel('Fidelity');

figure(4);
surf(dephase_exp,relax_exp,fidelity_std);
title('Fidelity std');
xlabel('log_{10}(dephase)');
ylabel('log_{10}(relax)');
zlabel('Fidelity');