%%insted the constraints ar*X-=Br + X*ac=Bc we have the dual constrain map(X) = am.
%%output-
%map - map(m) is all the linear eqution that we get by using ar*m and m*ac
%tmap = map'
%am = map(m)
%C - map(m) = C*vec(m).
function [map,tmap,am,C] = colrow_map(ar,ac,n,p,m,Br,Bc)
A = ar;
for i=2:n
    A = blkdiag(A,ar);
end

B = ac(1,1)*eye(n);
for j=2:n
    B = horzcat(B,ac(j,1)*eye(n));
end

for i=2:p
    C = ac(1,i)*eye(n);
    for j=2:n
        C = horzcat(C,ac(j,i)*eye(n));
    end
    
    B = vertcat(B,C);
end

C = vertcat(A,B);

am = vertcat(Br(:),Bc(:));
map = @(x)C*x(:);
tmap = @(x)C'*x;
end