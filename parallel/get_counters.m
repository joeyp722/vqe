function output = get_counters(count,count_array)
%Converts single counter to multiple counters
%Here count is the count number from the single counter, count_array the
%maximums for the multiple counter and output the count number for the
%mulitple counters.

product=1;
sum=0;
for i=1:length(count_array)
    counters(i)=(mod(count,count_array(i)*product)-sum)/product;
    sum=sum+counters(i)*product;
    product=product*count_array(i);
end

output=counters;