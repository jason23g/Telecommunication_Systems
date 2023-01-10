function X = bits_to_4_PAM(b,A)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i = 1 : 2: length(b) 
    
    if(b(i) == 0 && b(i+1) == 0)
        X((i + 1)/2) = 3*A;
    elseif(b(i) == 0 && b(i+1) == 1)
        X((i + 1)/2) = A;
    elseif(b(i) == 1 && b(i+1) == 1)
         X((i + 1)/2) = -A;
    elseif(b(i) == 1 && b(i+1) == 0)
        X((i + 1)/2) = -3*A;
  end
end

end

