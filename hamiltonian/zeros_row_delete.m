function output = zeros_row_delete(matrix)
%removes all the zero rows and colums
%input: matrix
%output: output

size=length(matrix);
row_delete=[];

for i = 1:size
    if sum(abs(matrix(i,:)))==0
        row_delete(end+1)=i;
    end
end

matrix(row_delete,:)=[];

matrix=matrix';

size=length(matrix);
row_delete=[];

for i = 1:size
    if sum(abs(matrix(i,:)))==0
        row_delete(end+1)=i;
    end
end

matrix(row_delete,:)=[];

matrix=matrix';

size=length(matrix);
nqbits_red=ceil(log(size)/log(2));

output_matrix=zeros(2^nqbits_red,2^nqbits_red);
output_matrix(1:size,1:size)=matrix;

output=output_matrix;



