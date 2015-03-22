

 function AX = Amap_Lu(X,G); 

 if isstruct(X)
    tmp = (G*X.U)*X.V';    
 else
    tmp = G*X; 
 end
 AX = tmp(:); 
