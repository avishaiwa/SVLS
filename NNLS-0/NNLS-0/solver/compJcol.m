%%***************************************************
%% find the number of nonzero elements (cumulatively) 
%% in each column of a sparse matrix
%%
%% NNLS, version 0: 
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun 
%%***************************************************

   function jcB = compJcol(JJ); 

   jcBtmp = [0; find(diff(JJ))];
   coltmp = JJ(jcBtmp+1); 
   cnt = 0;  
   for k = 1:length(jcBtmp) 
       if (k==1) 
          len = coltmp(k);
       else
          len = coltmp(k)-coltmp(k-1); 
       end
      jcB(cnt+[1:len]) = jcBtmp(k)*ones(1,len);  
      cnt = cnt + len; 
   end
   jcB = [jcB,length(JJ)]; 
%%***************************************************
