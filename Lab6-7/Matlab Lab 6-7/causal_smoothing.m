% Causal Smoothing Filter
function result =  causal_smoothing (a_signal,alpha , x2_k , s)

y4_k = zeros(length (x2_k ),1) ;

for i =  3 : length(x2_k)  ;
    
    y4_k(i) =  -s * alpha *  a_signal *  x2_k(i - 1) + (2*a_signal)* y4_k(i-1) - (a_signal^2) * y4_k(i-2) ;
    
end

result = y4_k ;

end