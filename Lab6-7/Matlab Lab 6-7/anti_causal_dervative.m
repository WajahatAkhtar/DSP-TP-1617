%Anti - Causal Sharpening Filter
function result = anti_causal_dervative (a_signal,alpha , x_k)

y_k = zeros(length (x_k ),1) ;
R = length(x_k)-2 : -1 : 1 ;

for i =  R
    
    y_k(i) =  a_signal   *   ( alpha + 1 ) *  x_k(i+1) - (a_signal^2) * x_k(i+2) + (2*a_signal)* y_k(i+1) - (a_signal^2) * y_k(i+2) ;
    
end

result = y_k;

end