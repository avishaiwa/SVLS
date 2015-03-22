

  function Aty = ATmap_Lu(y,G); 
  
  
  [nr,nc] = size(G); 
  p = length(y)/nr; 
  Y = reshape(y,nr,p); 
  Aty = G'*Y; 
