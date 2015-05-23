% Draw at random row and column measurements for SVLS algorithm 
function [ Br,Bc,ar,ac ] =create_measurements( n,k,am )
if length(n) == 1,
    n1 = n(1); n2 = n1;
elseif length(n) == 2,
    n1 = n(1); n2 = n(2);
end

col_samp=randsample(n2,k);
row_samp=randsample(n1,k);

k_col=1;
k_row = 1;
%Create Ac and Bc
for i=1:n2
    if (length(find(am(:,i)~=0))==n1)
        col_samp(k_col) = i;
        k_col= k_col +1;
    end
end
%Create Ar and Br
for i=1:n1
      if (length(find(am(i,:))~=0)==n2)
        row_samp (k_row) = i;
         k_row= k_row +1;
      end
end

ar = zeros(k,n1);
ac = zeros(n2,k);
for i=1:k
        ar(i,row_samp(i))=1;
end
for i=1:k
        ac(col_samp(i),i)=1;
end

Br = ar*am;
Bc = am*ac;
end
