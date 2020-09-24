function RhoOut = EntanglementGateSelect(RhoIn,param)
%Selects entanglement gate sequence based on the implementation type.
%input: RhoIn, param struct
%these are:
%input density matrix, parameter for implementation type
%output: RhoOut
%this is: output density matrix

Rho=RhoIn;

switch param.select
    case 'Ion'
        Rho=IonEntanglementGate(Rho,param);
    case 'Majorana'
        Rho=MajoranaEntanglementGate(Rho,param);
     case 'Silicon'
        Rho=SiliconEntanglementGate(Rho,param);
    case 'Transmon'
        Rho=TransmonEntanglementGate(Rho,param);
    case 'General'
        Rho=GeneralEntanglementGate(Rho,param);

end

RhoOut = Rho;