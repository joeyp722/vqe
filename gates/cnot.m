function output = cnot(control,target,nqbits)
%Generates CNOT matrix for specified control and target qubits.
%input: control, target, nqbits
%these are:
%the control qubit, target qubit, total number of qubits
%output: output
%this is:
%the cnot matrix with properties stated above

size=2^nqbits;
cnot=zeros(size,size);

for i=1:size
    Q=i-1;
    bin_vec=decimalToBinaryVector(Q,nqbits);

    if(bin_vec(control)==1)
        bin_vec(target)=~bin_vec(target);
        Q=binaryVectorToDecimal(bin_vec);
        cnot(i,Q+1)=1;
    else
        cnot(i,Q+1)=1;
    end      
end

output=cnot;