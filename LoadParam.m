function param = LoadParam(select)
%Loads parameters for the type of physical implementation that is
%considered.
%input: select
%this is:
%name of implementation
%output: param
%these are: The parameters related to the physical cosntants of the
%implementation.

switch select
    case 'Ion'
        param.rabi=3.2*10^9; %Hz
        param.detuning=0.4*10^9; %Hz
        
        phase_whitenoise_carrier_ratio=-152; %dBc
        amplitude_whitenoise_carrier_ratio=-140; %dBc
        
        param.dephase=1/(50*10^9); %Hz
        param.amplitude=param.rabi^2*10^(amplitude_whitenoise_carrier_ratio/10)/8;
        param.phase=param.rabi^2*10^(phase_whitenoise_carrier_ratio/10)/8;
     
        param.entanglement_time=0;
        
    case 'Majorana'
        param.quasiparticle.p=10^-3;
        param.quasiparticle.p_qp=1.6*10^-5;
        param.quasiparticle.p_cor_even=3.2*10^-5;
        param.quasiparticle.p_cor_odd=8*10^-6;
        param.quasiparticle.p_pair=6.4*10^-5;
        
    case 'Silicon'
        param.rabi=13*10^9; %Hz
        
        param.relax=1/(0.5*10^-6); %Hz
        param.dephase=1/(37*10^-6); %Hz
                
        param.entanglement_time=3.3*10^-6; %seconds
        
    case 'Transmon'
        param.rabi=5.3*10^9; %Hz
        
        param.relax=1/(30*10^6); %Hz
        param.dephase=1/(30*10^6); %Hz
                
        param.entanglement_time=450*10^-9; %seconds
        
	case 'General'
        param.rabi=1;
        
        param.relax=1;
        param.dephase=1;
                
        param.entanglement_time=1;
end