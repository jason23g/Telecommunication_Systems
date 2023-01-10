function N = compare_symbols(X,Y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

counter = 0;

for i = 1:length(X)
    counter1 = 0;
    if(X(i) ~= Y(i))
        counter1 = counter + 1;
        counter = counter1;
    end
    
end
N = counter;
end

