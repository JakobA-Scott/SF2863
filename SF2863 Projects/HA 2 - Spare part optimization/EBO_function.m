function [EBO] = EBO_function(S_n, n)

Lambda = [50 40 45 50 25 48 60 35 15]/1000; %Arrival rates of malfunctioning LRU's
T = [6 8 14 25 12 18 33 8 12];  %Repair time for malfunctioning LRU
c = [12 14 21 20 11 45 75 30 22];   %Purchase cost per LRU


EBO = Lambda(n)*T(n);
state = 0;
s = state;
%We are sending in the vector S_n
state_check = S_n>0;

if state_check == true   

   while s < S_n
       
        EBO = EBO - (1 - poisscdf(s, Lambda(n)*T(n)));
        s = s+1;
        
    end

end

end
