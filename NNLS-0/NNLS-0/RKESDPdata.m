%%*************************************************************************
%% generate SDP data corresponding to 
%%
%% min_{X psd} sum_{ij\in nnz(W)} (<Aij,X>-dij)^2 + lam*Tr(X)
%%
%% input: DD = (npts)x(npts) dis-similar matrix 
%% 
%%*************************************************************************

    function  [blk,At,bb,numeq] = RKESDPdata(DD,WW); 
  
    if (nargin < 2); WW = spones(DD); end
    [ii,jj,vv] = find(WW); 
    if (min(vv) < 0); error('nonzero element of WW must be positive'); end
    npts = size(DD,2);
    if (size(DD,1) ~= size(DD,2)); error(' DD must be square'); end
    DD = triu(DD,1);
    DD = DD.*spones(WW);  
    NZ = nnz(DD);
%%     
    sqrt2 = sqrt(2); 
    bb = zeros(NZ,1); 
    II = zeros(3*NZ,1);
    JJ = zeros(3*NZ,1);
    VV = zeros(3*NZ,1); 
    count = 0; 
    mm = 0; 
    for k = 1:npts
       for j = 1:k-1
          rr = DD(j,k);
          if (rr)
             mm = mm+1;  
             bb(mm) = rr;
	     II(count+[1:3]) = [j*(j+1)/2;(k-1)*k/2+j;k*(k+1)/2];
             JJ(count+[1:3]) = mm*ones(3,1);
             VV(count+[1:3]) = [1;-sqrt2;1]; 
             count = count + 3; 
          end
       end
    end
%%
%% set e'*Y*e = 0 to center points around origin
%%
    blk{1,1} = 's'; blk{1,2} = npts;    
    ee = sqrt(1/npts)*ones(npts,1); 
    Alast = svec(blk,ee*ee',1); 
    At{1} = [spconvert([II,JJ,VV;npts*(npts+1)/2,mm,0]),Alast];
    bb(mm+1) = 0; 
    numeq = 1; 
%%************************************************************************
