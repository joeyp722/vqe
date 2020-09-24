function RhoOut = MultiRotationGateSelect(RhoIn,param,angles)
%Selects rotational gate sequence based to the implementation type
%input: RhoIn, param, angles
%these are:
%input density matrix, selection of implemention type, rotational angles of
%the gates (this is a matrix)
%output: RhoOut
%this is: output density matrix

Rho=RhoIn;

switch param.select
    case 'Ion'
        Rho=MultiIonRotationGate(Rho,param,angles);
	case 'Majorana'
        Rho=MultiMajoranaRotationGate(Rho,param,angles);
	case 'Silicon'
        Rho=MultiSiliconRotationGate(Rho,param,angles);
    case 'Transmon'
        Rho=MultiTransmonRotationGate(Rho,param,angles);
    case 'General'
        Rho=MultiGeneralRotationGate(Rho,param,angles);
end

RhoOut=Rho;