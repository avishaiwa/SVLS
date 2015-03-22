%%**********************************************************************
%% find the vector obtained from the matrix
%% X = U*V' by extracting the entries with indices specified 
%% in (II,Jcol), which is in compressed sparse column format. 
%%
%% NNLS, version 0: 
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun 
%%**********************************************************************

   function xx = ProjOmega(U,V,II,Jcol); 

   if (exist('mexProjOmega') == 3)
      xx = mexProjOmega(U,V,II,Jcol); 
   else
      nn = length(II); 
      xx = zeros(nn,1); 
      len = length(Jcol)-1; 
      cnt  = 0; 
      if issparse(U); U = full(U); V = full(V); end
      for k = 1:len
         VTk = V(k,:)';
         Xk = U*VTk; 
         idx = [Jcol(k)+1 : Jcol(k+1)]; 
         xx(idx) = Xk(II(idx));       
         cnt = cnt + length(idx); 
      end
   end
%%**********************************************************************

