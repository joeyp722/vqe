%example script calculate the potentials of all states for molecular hydrogen considering the ion
%implementation without errors, here the ssvqe weighted algorithm is used.

clear all;

%select implemtation type
select='Ion';
param = LoadParam(select);
param.select=select;

%set iterations for repition of the ssvqe subroutine
iterations=5;

%set ssvqe subroute parameters
param.steps=400; %iterations until solution is found
param.weight_multiplier=10; %weight multiplier for cost function
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
param.error=false;

%ditance parameters
param.distance_interval=0.05;
param.min_distance=0.10;
param.max_distance=2.30;

distance_interval=param.distance_interval;
min_distance=param.min_distance;
max_distance=param.max_distance;

max_distance_index=ceil((max_distance-min_distance+distance_interval)/distance_interval);

name=param.name;

for distance_index=1:max_distance_index
    
    %set distance
    param.distance=(distance_index-1)*distance_interval+min_distance;
    
    %load matrices
    load(strcat('molecules/',name,'/one_electron_matrix',num2str(distance_index),'.mat'));
    one_electron_matrix=matrix;
    load(strcat('molecules/',name,'/two_electron_matrix',num2str(distance_index),'.mat'));
    two_electron_matrix=matrix;

    %construct hamiltonian
    param.hamiltonian=get_hamiltonian(one_electron_matrix,two_electron_matrix,param);
    nstates=length(param.hamiltonian); %get the number of states
    
    
    for i=1:iterations
        %ssvqe subroutine
        out=ssvqe_weight(param);
        output{i,distance_index}=out;
    end

end

%calculating the mean and std of the potentials, fidelities and pseudo potentials for multiple states
%based on multiple iterations
for j=1:nstates
    for distance_index=1:max_distance_index
        potential_dataset=[];
        fidelity_dataset=[];
        pseudo_potential_dataset=[];
        for i=1:iterations
            potential_dataset(i)=output{i,distance_index}.Potential{j};
            fidelity_dataset(i)=output{i,distance_index}.Fidelity{j};
            pseudo_potential_dataset(i)=output{i,distance_index}.pseudo_Potential{j};
        end
        data{distance_index}.potential.mean{j}=mean(potential_dataset);
        data{distance_index}.potential.std{j}=std(potential_dataset);
        data{distance_index}.fidelity.mean{j}=mean(fidelity_dataset);
        data{distance_index}.fidelity.std{j}=std(fidelity_dataset);
        data{distance_index}.pseudo_potential.mean{j}=mean(pseudo_potential_dataset);
        data{distance_index}.pseudo_potential.std{j}=std(pseudo_potential_dataset);
    end
end

%defining plotting axes
distance_axis=[];
potential_curve_mean_axis=[];
potential_curve_std_axis=[];
fidelity_curve_mean_axis=[];
fidelity_curve_std_axis=[];
pseudo_potential_curve_mean_axis=[];
pseudo_potential_curve_std_axis=[];

%Putting the data in the plotting axes
for j=1:nstates
    for distance_index=1:max_distance_index
        distance_axis(distance_index)=(distance_index-1)*distance_interval+min_distance;
        potential_curve_mean_axis(distance_index,j)=data{distance_index}.potential.mean{j};
        potential_curve_std_axis(distance_index,j)=data{distance_index}.potential.std{j};
        fidelity_curve_mean_axis(distance_index,j)=data{distance_index}.fidelity.mean{j};
        fidelity_curve_std_axis(distance_index,j)=data{distance_index}.fidelity.std{j};
        pseudo_potential_curve_mean_axis(distance_index,j)=data{distance_index}.pseudo_potential.mean{j};
        pseudo_potential_curve_std_axis(distance_index,j)=data{distance_index}.pseudo_potential.std{j};
    end

    %plotting results
    
    legend_states{j}=strcat(string(j-1),"-th excited state");
    
    hold on
    
    figure(1);
    plot(distance_axis,potential_curve_mean_axis(:,j));
    title('Potential mean');
    xlabel('distance Angstroms');
    ylabel('Energy [Hartrees]');
    legend(legend_states);

    hold off
    hold on
    
    figure(2);
    plot(distance_axis,potential_curve_std_axis(:,j));
    title('Potential std');
    xlabel('distance Angstroms');
    ylabel('Energy [Hartrees]');
    legend(legend_states);
    
    hold off
    hold on
    
    figure(3);
    plot(distance_axis,fidelity_curve_mean_axis(:,j));
    title('Fidelity mean');
    xlabel('distance Angstroms');
    ylabel('Fidelity');
    legend(legend_states);
    
    hold off
    hold on
    
    figure(4);
    plot(distance_axis,fidelity_curve_std_axis(:,j));
    title('Fidelity std');
    xlabel('distance Angstroms');
    ylabel('Fidelity');
    legend(legend_states);
    
    hold off
    hold on
    
    figure(5);
    plot(distance_axis,pseudo_potential_curve_mean_axis(:,j));
    title('pseudo Potential mean');
    xlabel('distance Angstroms');
    ylabel('Energy [Hartrees]');
    legend(legend_states);

    hold off
    hold on
    
    figure(6);
    plot(distance_axis,pseudo_potential_curve_std_axis(:,j));
    title('pseudo Potential std');
    xlabel('distance Angstroms');
    ylabel('Energy [Hartrees]');
    legend(legend_states);
    
    hold off 
    hold on
    
    figure(7);
    plot(distance_axis,potential_curve_mean_axis(:,j)-pseudo_potential_curve_mean_axis(:,j));
    title('difference Potential mean');
    xlabel('distance Angstroms');
    ylabel('Energy [Hartrees]');
    legend(legend_states);

    hold off
    hold on
    
    figure(8);
    plot(distance_axis,potential_curve_std_axis(:,j)-pseudo_potential_curve_std_axis(:,j));
    title('difference Potential std');
    xlabel('distance Angstroms');
    ylabel('Energy [Hartrees]');
    legend(legend_states);
    
    hold off
    
end