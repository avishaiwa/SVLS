%%**********************************************************************
%% find X.*mask where X = U*V'; 
%% mask = spconvert([II,JJ,ones(length(II),1);nr nc 1]); 
%%
%% NNLS, version 0: 
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun 
%%**********************************************************************

   function xx = ProjOmega(U,V,II,JJ); 

   nn = length(II); 
   xx = zeros(nn,1); 
   Jcol = [0; find(diff(JJ)); nn];  
   len = length(Jcol)-1; 
   cnt  = 0; 
   if issparse(U); U = full(U); V = full(V); end
   for k = 1:len
      j  = JJ(Jcol(k)+1);
      Vj = V(j,:); 
      Xj = U*Vj'; 
      idx = [Jcol(k)+1 : Jcol(k+1)]; 
      xx(idx) = Xj(II(idx));       
      cnt = cnt + length(idx); 
   end
%%**********************************************************************

