function Ionization_Energies = get_Ionization_Energies(name)
%Gives ionization energy for molecules stated below
%input: name, distance
%these are: name of molecule, distance between atoms
%output: Ionization_Energies
%this is: the ionization energy

switch name
    case 'H2'
        Ionization_Energies=0;
    case 'LiH'
        Ionization_Energies=-(75.6400964+122.4543581)/27.2;
    case 'LiH_core'
        Ionization_Energies=0;
    case 'LiH_active3'
        Ionization_Energies=-(75.6400964+122.4543581)/27.2;
    case 'BeH2'
        Ionization_Energies=-(153.896203+217.7185843)/27.2;
    case 'BeH2_active4'
        Ionization_Energies=-(153.896203+217.7185843)/27.2;
end