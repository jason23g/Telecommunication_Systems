function X = bits_to_2PAM(b)
	
	for k = 1:length(b)
		if (b(k) == 0)
			X(k) = 1;
		else
			X(k) = -1;
		end
	end
end