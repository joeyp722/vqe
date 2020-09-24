function RhoOut = LindbladSelect(RhoIn,param,gate)
%Executes relevant Kraus transformations for the type of implementation and
%gate

%input: RhoIn, param struct, gate
%these are: input density matrix, parameters defining the implementation
%and parameters relevant for the Kraus transformations, selection of
%rotational or entanglement gate

%output: RhoOut
%this is: output density matrix

Rho = RhoIn;

switch param.select
    case 'Ion'
        switch gate
            case 'rotational'
                dephase=param.dephase;
                amplitude=param.amplitude;
                phase=param.phase;
        
                rotational_time=param.rotational_time;
                
                Rho=LindbladDephase(Rho,rotational_time,dephase);
                Rho=LindbladAmplitude(Rho,rotational_time,amplitude);
                Rho=LindbladPhase(Rho,rotational_time,phase);
                
            case 'entanglement'
                dephase=param.dephase;
                amplitude=param.amplitude;
                phase=param.phase;
        
                entanglement_time=param.entanglement_time;
                
                Rho=LindbladDephase(Rho,entanglement_time,dephase);
                Rho=LindbladAmplitude(Rho,entanglement_time,amplitude);
                Rho=LindbladPhase(Rho,entanglement_time,phase);
        end

    case 'Majorana'
        switch gate
            case 'rotational'
                p=param.quasiparticle.p;
                p_qp=param.quasiparticle.p_qp;
                p_cor_even=param.quasiparticle.p_cor_even;
                p_cor_odd=param.quasiparticle.p_cor_odd;
                p_pair=param.quasiparticle.p_pair;
        
                Rho=LindbladQuasiParticle(Rho,p,p_qp,p_cor_even,p_cor_odd,p_pair);
                
            case 'entanglement'
                p=param.quasiparticle.p;
                p_qp=param.quasiparticle.p_qp;
                p_cor_even=param.quasiparticle.p_cor_even;
                p_cor_odd=param.quasiparticle.p_cor_odd;
                p_pair=param.quasiparticle.p_pair;
        
                Rho=LindbladQuasiParticle(Rho,p,p_qp,p_cor_even,p_cor_odd,p_pair);
        end
        
    case 'Silicon'
        switch gate
            case 'rotational'
                relax=param.relax;
                dephase=param.dephase;
                
                rotational_time=param.rotational_time;
        
                Rho=LindbladRelax(Rho,rotational_time,relax);
                Rho=LindbladDephase(Rho,rotational_time,dephase);
                
            case 'entanglement'
                relax=param.relax;
                dephase=param.dephase;
                
                entanglement_time=param.entanglement_time;
        
                Rho=LindbladRelax(Rho,entanglement_time,relax);
                Rho=LindbladDephase(Rho,entanglement_time,dephase);
        end
        
	case 'Transmon'
        switch gate
            case 'rotational'
                relax=param.relax;
                dephase=param.dephase;
                
                rotational_time=param.rotational_time;
        
                Rho=LindbladRelax(Rho,rotational_time,relax);
                Rho=LindbladDephase(Rho,rotational_time,dephase);
                
            case 'entanglement'
                relax=param.relax;
                dephase=param.dephase;
                
                entanglement_time=param.entanglement_time;
        
                Rho=LindbladRelax(Rho,entanglement_time,relax);
                Rho=LindbladDephase(Rho,entanglement_time,dephase);
        end
        
	case 'General'
        switch gate
            case 'rotational'
                relax=param.relax;
                dephase=param.dephase;
                
                rotational_time=param.rotational_time;
        
                Rho=LindbladRelax(Rho,rotational_time,relax);
                Rho=LindbladDephase(Rho,rotational_time,dephase);
                
            case 'entanglement'
                relax=param.relax;
                dephase=param.dephase;
                
                entanglement_time=param.entanglement_time;
        
                Rho=LindbladRelax(Rho,entanglement_time,relax);
                Rho=LindbladDephase(Rho,entanglement_time,dephase);
        end
end

RhoOut=Rho;