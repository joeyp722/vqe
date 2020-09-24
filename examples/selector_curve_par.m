%example script calculate the potential of the ground state for molecular hydrogen considering the ion
%implementation. With parallel processing.

clear all;

%select implemtation type
select='Ion';
param = LoadParam(select);
param.select=select;

%set iterations for repition of the vqe subroutine
iterations=5;

%set vqe subroute parameters
param.steps=5000; %iterations until solution is found
param.depth=3; %circuit depth

%SPSA parameters
param.alpha=0.602;
param.gamma=0.101;
param.c=0.1;

%parameters for construct hamiltonian
param.name='H2'; %molecule name
param.transformation_type='jw red'; %hamiltonian construction method
param.electron_number=2; %number of valence electrons
param.spin_number=0; %total valence spin number

%error flag
param.error=true;

%ditance parameters
param.distance_interval=0.05;
param.min_distance=0.10;
param.max_distance=5.00;
distance_interval=param.distance_interval;
min_distance=param.min_distance;
max_distance=param.max_distance;

max_distance_index=ceil((max_distance-min_distance+distance_interval)/distance_interval);

name=param.name;

%Prepare electron matrices for each parallel process.
for distance_index=1:max_distance_index
    load(strcat('molecules/',name,'/one_electron_matrix',num2str(distance_index),'.mat'));
    one_electron_matrix{distance_index}=matrix;
    load(strcat('molecules/',name,'/two_electron_matrix',num2str(distance_index),'.mat'));
    two_electron_matrix{distance_index}=matrix;
    
    %Constructing param for each distance
    param_dis{distance_index}=param;
    
    %set distance
    param_dis{distance_index}.distance=(distance_index-1)*distance_interval+min_distance;
    
    %construct hamiltonian
    param_dis{distance_index}.hamiltonian=get_hamiltonian(one_electron_matrix{distance_index},two_electron_matrix{distance_index},param);
    
    %Display progress
    disp(distance_index);
end

%set counter variables
count_array=[max_distance_index,iterations];
count_total=prod(count_array);

parfor_progress(count_total);
parfor j=1:count_total
    
    %counter
    counters=get_counters(j-1,count_array);
    i_count=counters(2);
    distance_count=counters(1);
   
    %counter variables conversion
    distance_index=distance_count+1;
    i=i_count+1;
    
    %Constructing param for each parallel process
    param_par(j)=param_dis{distance_index};
    
    param_par(j).error=false;
        
    %vqe subroutine without errors
    output_par(j).noerror=vqe(param_par(j));
        
    param(j).error=true;
        
    %vqe subroutine with errors
    output_par(j).error=vqe(param_par(j));

    parfor_progress;
end
parfor_progress(0);

%Unpack parallel processing results
for j=1:count_total
    
    %counter
    counters=get_counters(j-1,count_array);
    i_count=counters(2);
    distance_count=counters(1);
   
    %counter variables conversion
    distance_index=distance_count+1;
    i=i_count+1;
    
    output{i}.noerror{distance_index}=output_par(j).noerror;
    output{i}.error{distance_index}=output_par(j).error;
    
end

%calculating the mean and std of the potential based on multiple iterations without errors
for distance_index=1:max_distance_index
    dataset=[];
    for i=1:iterations
        dataset(i)=output{i}.noerror{distance_index}.Potential;
    end
    output{1}.noerror_potential{distance_index}.mean=mean(dataset);
    output{1}.noerror_potential{distance_index}.std=std(dataset);
end

%calculating the mean and std of the potential based on multiple iterations with errors
for distance_index=1:max_distance_index
    dataset=[];
    for i=1:iterations
        dataset(i)=output{i}.error{distance_index}.Potential;
    end
    output{1}.error_potential{distance_index}.mean=mean(dataset);
    output{1}.error_potential{distance_index}.std=std(dataset);
end

%calculating the mean and std of the fidelity based on multiple iterations without errors
for distance_index=1:max_distance_index
    dataset=[];
    for i=1:iterations
        dataset(i)=output{i}.noerror{distance_index}.Fidelity;
    end
    output{1}.noerror_fidelity{distance_index}.mean=mean(dataset);
    output{1}.noerror_fidelity{distance_index}.std=std(dataset);
end

%calculating the mean and std of the fidelity based on multiple iterations without errors
for distance_index=1:max_distance_index
    dataset=[];
    for i=1:iterations
        dataset(i)=output{i}.error{distance_index}.Fidelity;
    end
    output{1}.error_fidelity{distance_index}.mean=mean(dataset);
    output{1}.error_fidelity{distance_index}.std=std(dataset);
end

%calculating the mean and std of the pseudo potential based on multiple iterations without errors
for distance_index=1:max_distance_index
    dataset=[];
    for i=1:iterations
        dataset(i)=output{i}.noerror{distance_index}.pseudo_Potential;
    end
    output{1}.pseudo_potential{distance_index}.mean=mean(dataset);
    output{1}.pseudo_potential{distance_index}.std=std(dataset);
end

%constructing plotting data
distance_axis=[];
potential_noerror_mean_axis=[];
potential_error_mean_axis=[];
potential_noerror_std_axis=[];
potential_error_std_axis=[];
fidelity_noerror_mean_axis=[];
fidelity_error_mean_axis=[];
fidelity_noerror_std_axis=[];
fidelity_error_std_axis=[];
pseudo_potential_mean_axis=[];
pseudo_potential_std_axis=[];

for distance_index=1:max_distance_index
distance_axis(distance_index)=(distance_index-1)*distance_interval+min_distance;
potential_noerror_mean_axis(distance_index)=output{1}.noerror_potential{distance_index}.mean;
potential_error_mean_axis(distance_index)=output{1}.error_potential{distance_index}.mean;
potential_noerror_std_axis(distance_index)=output{1}.noerror_potential{distance_index}.std;
potential_error_std_axis(distance_index)=output{1}.error_potential{distance_index}.std;
fidelity_noerror_mean_axis(distance_index)=output{1}.noerror_fidelity{distance_index}.mean;
fidelity_error_mean_axis(distance_index)=output{1}.error_fidelity{distance_index}.mean;
fidelity_noerror_std_axis(distance_index)=output{1}.noerror_fidelity{distance_index}.std;
fidelity_error_std_axis(distance_index)=output{1}.error_fidelity{distance_index}.std;
pseudo_potential_mean_axis(distance_index)=output{1}.pseudo_potential{distance_index}.mean;
pseudo_potential_std_axis(distance_index)=output{1}.pseudo_potential{distance_index}.std;
end

difference_potential_noerror_mean_axis=potential_noerror_mean_axis-pseudo_potential_mean_axis;
difference_potential_noerror_std_axis=potential_noerror_std_axis-pseudo_potential_std_axis;
difference_potential_error_mean_axis=potential_error_mean_axis-pseudo_potential_mean_axis;
difference_potential_error_std_axis=potential_error_std_axis-pseudo_potential_std_axis;

%plotting data
figure(1);
plot(distance_axis,pseudo_potential_mean_axis,distance_axis,potential_noerror_mean_axis,distance_axis,potential_error_mean_axis);
title('Potential mean');
xlabel('distance Angstroms');
ylabel('Energy [Hartrees]');
legend('pseudo','no error','error');

figure(2);
plot(distance_axis,pseudo_potential_std_axis,distance_axis,potential_noerror_std_axis,distance_axis,potential_error_std_axis);
title('Potential std');
xlabel('distance Angstroms');
ylabel('Energy [Hartrees]');
legend('pseudo','no error','error');

figure(3);
plot(distance_axis,difference_potential_noerror_mean_axis,distance_axis,difference_potential_error_mean_axis);
title('difference Potential mean');
xlabel('distance Angstroms');
ylabel('Energy [Hartrees]');
legend('noerror','error');

figure(4);
plot(distance_axis,difference_potential_noerror_std_axis,distance_axis,difference_potential_error_std_axis);
title('difference Potential std');
xlabel('distance Angstroms');
ylabel('Energy [Hartrees]');
legend('noerror','error');

figure(5);
plot(distance_axis,fidelity_noerror_mean_axis,distance_axis,fidelity_error_mean_axis);
title('Fidelity mean');
xlabel('distance Angstroms');
ylabel('Fidelity');
legend('no error','error');

figure(6);
plot(distance_axis,fidelity_noerror_std_axis,distance_axis,fidelity_error_std_axis);
title('Fidelity std');
xlabel('distance Angstroms');
ylabel('Fidelity');
legend('no error','error');