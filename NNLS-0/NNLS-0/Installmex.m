%%***********************************************************************
%% compile mex files
%%
%% NNLS, version 0: 
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun 
%%***********************************************************************

   function Installmex

   warning off
   computer_model = computer;
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);
%%
   if strcmp(computer_model,'PCWIN')
      str1 = ['  ''',matlabroot,'\extern\lib\win32\lcc\libmwlapack.lib''  ']; 
      str2 = ['  ''',matlabroot,'\extern\lib\win32\lcc\libmwblas.lib''  '];
      libstr = [str1,str2];     
   elseif strcmp(computer_model,'PCWIN64')
      str1 = ['  ''',matlabroot,'\extern\lib\win64\microsoft\libmwlapack.lib''  '];
      str2 = ['  ''',matlabroot,'\extern\lib\win64\microsoft\libmwblas.lib''  '];
      libstr = [str1,str2];  
   else
      libstr = ' -lmwlapack -lmwblas  '; 
   end
   mexcmd = 'mex -O  -largeArrayDims  -output ';    

   cd PROPACKmod
   cmd([mexcmd, 'reorth mexreorth.c','   ',libstr]);
   cd ..
  
   cd solver
   if (ispc)
      cmd([mexcmd, 'mexeig mexeigwin.c','  ',libstr]);
      cmd([mexcmd, 'mexsvd mexsvdwin.c','  ',libstr]);
   else
      cmd([mexcmd, 'mexeig mexeig.c','  ',libstr]);  
      cmd([mexcmd, 'mexsvd mexsvd.c','  ',libstr]);  
   end
   fname{1} = 'mexsvec'; 
   fname{2} = 'mexsmat'; 
   fname{3} = 'mexProjOmega';
   fname{4} = 'mexspconvert'; 
   for k = 1:length(fname)
      cmd([mexcmd,'  ',fname{k},'  ',fname{k},'.c','  ',libstr]);    
   end
   cd ..
%%***********************************************
   function cmd(s) 
   
   fprintf(' %s\n',s); 
   eval(s); 
%%***********************************************
