%Anti - Causal Smoothing Filter 
function result =  anticausal_smoothing (a_signal,alpha , x2_k ,s)

y3_k = zeros(length (x2_k ),1) ;
R = length(x2_k)-2 : -1 : 1 ;

for i =  R
    
    y3_k(i) = s * alpha *  a_signal   *   x2_k(i+1)  + (2*a_signal)* y3_k(i+1) - (a_signal^2) * y3_k(i+2) ;
    
end

result = y3_k ;

end