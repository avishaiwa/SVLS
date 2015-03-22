%%**************************************************
%% compute A*z for 
%% A = beta1*X - beta2*Xold - beta3*Grad;
%% 
%% Az = matvec(z,X,Xold,Grad,beta1,beta2,tau); 
%%
%% NNLS, version 0: 
%% Copyright (c) 2009 by
%% Kim-Chuan Toh and Sangwoon Yun 
%%**************************************************

  function Az = Amatvec(z,param); 

  global Xh Xhold Grad 

  beta1 = param.beta1; 
  beta2 = param.beta2; 
  beta3 = param.beta3; 

  z1 = Xh.U*(Xh.V'*z); 
  z2 = Xhold.U*(Xhold.V'*z); 
  Az = beta1*z1 -beta2*z2 -beta3*(Grad*z); 
%%**************************************************
