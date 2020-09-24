function Repulsion = get_Repulsion(name,distance)
%Gives repulsion energy for molecules stated below
%input: name, distance
%these are: name of molecule, distance between atoms
%output: Repulsion
%this is: the repulsion energy

switch name
    case 'H2'
        Repulsion=(0.03674932*1.6021766E-19)/(4*pi*8.85418782E-12)/(distance*10^-10);
    case 'LiH'
        Z=1.279;
        Repulsion=(Z*0.03674932*1.6021766E-19)/(4*pi*8.85418782E-12)/(distance*10^-10);
    case 'LiH_core'
        Z=3;
        Repulsion=(Z*0.03674932*1.6021766E-19)/(4*pi*8.85418782E-12)/(distance*10^-10);
    case 'LiH_active3'
        Z=1.279;
        Repulsion=(Z*0.03674932*1.6021766E-19)/(4*pi*8.85418782E-12)/(distance*10^-10);
    case 'BeH2'
        Z=1.912;
        Repulsion=2*(Z*0.03674932*1.6021766E-19)/(4*pi*8.85418782E-12)/(distance*10^-10)+(0.03674932*1.6021766E-19)/(4*pi*8.85418782E-12)/(2*distance*10^-10);
    case 'BeH2_active4'
        Z=1.912;
        Repulsion=2*(Z*0.03674932*1.6021766E-19)/(4*pi*8.85418782E-12)/(distance*10^-10)+(0.03674932*1.6021766E-19)/(4*pi*8.85418782E-12)/(2*distance*10^-10);
end