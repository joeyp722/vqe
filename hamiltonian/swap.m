function output = swap(H)
%swaps ground and excited state in a Hamiltonian
%input: H
%this is: input Hamiltonian
%output: output
%this is: output Hamiltonian

size=length(H);
swap_matrix=zeros(size,size);
for i=1:size
    swap_matrix(i,size-i+1)=1;
end

output=swap_matrix*H*swap_matrix;