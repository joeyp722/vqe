function output=get_hamiltonian(one_electron_matrix,two_electron_matrix,param)
%Constructs the hamiltonain in Pauli basis based on the one and two
%electron matrices, number of electrons, spin number and transformation
%type. These types are the Jordan-Wigner transformation 'jw' and the JW
%with qubit reduction 'jw red'

%input: one_electron_matrix, two_electron_matrix, param
%these are: total one electron matrix, total two electron matrix, struct
%where the additional parameters above are defined
%output: output
%this is: the construct Hamiltonian matrix

transformation_type=param.transformation_type;
electron_number=param.electron_number;
spin_number=param.spin_number;

size_one=length(one_electron_matrix);
size_two=length(two_electron_matrix);

one_electron_hamiltonian=0;
two_electron_hamiltonian=0;

switch transformation_type
    case 'jw'
        for i=1:size_one
            for j=1:size_one
                one_electron_hamiltonian=one_electron_hamiltonian+one_electron_matrix(i,j)*jordanwigner(i,size_one,-1)*jordanwigner(j,size_one,1);
            end
        end
        
        for i=1:size_two
            for j=1:size_two
                for k=1:size_two
                    for l=1:size_two
                        two_electron_hamiltonian=two_electron_hamiltonian+(1/2)*two_electron_matrix(i,j,k,l)*jordanwigner(i,size_two,-1)*jordanwigner(j,size_two,-1)*jordanwigner(k,size_two,1)*jordanwigner(l,size_two,1);
                    end
                end
            end
        end
        hamiltonian=one_electron_hamiltonian+two_electron_hamiltonian;
        output=hamiltonian;
        
    case 'jw red'
        for i=1:size_one
            for j=1:size_one
                one_electron_hamiltonian=one_electron_hamiltonian+one_electron_matrix(i,j)*jordanwigner(i,size_one,-1)*jordanwigner(j,size_one,1);
            end
        end
        
        for i=1:size_two
            for j=1:size_two
                for k=1:size_two
                    for l=1:size_two
                        two_electron_hamiltonian=two_electron_hamiltonian+(1/2)*two_electron_matrix(i,j,k,l)*jordanwigner(i,size_two,-1)*jordanwigner(j,size_two,-1)*jordanwigner(k,size_two,1)*jordanwigner(l,size_two,1);
                    end
                end
            end
        end
        hamiltonian=one_electron_hamiltonian+two_electron_hamiltonian;
        
        nqbits=log(length(hamiltonian))/log(2);
        hamiltonian_projected=transpose(projection(electron_number,nqbits))*hamiltonian*projection(electron_number,nqbits);
        hamiltonian_spinprojected=transpose(spinprojection(spin_number,nqbits))*hamiltonian_projected*spinprojection(spin_number,nqbits);
        hamiltonian_red=zeros_row_delete(hamiltonian_spinprojected);
        
        hamiltonian=hamiltonian_red;
        
        output=hamiltonian;
        
 
    otherwise
        output=0;
end
