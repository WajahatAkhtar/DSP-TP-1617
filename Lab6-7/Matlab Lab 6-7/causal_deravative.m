% Causal Sharpening Filter  

function result =  causal_deravative (a_signal,alpha , x_k)
y2_k = zeros(length (x_k ),1) ;

for i = 3 : length(x_k)
    
    y2_k(i) =  x_k(i) + a_signal   *   ( alpha - 1 ) *  x_k(i - 1) + (2*a_signal)* y2_k(i-1) - (a_signal^2) * y2_k(i-2) ;
    
end

result = y2_k ;

end
