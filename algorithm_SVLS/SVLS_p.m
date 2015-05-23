function [X] = SVLS_p( Br,Bc,Ar,Ac,r, qr_flag)

if(~exist('qr_flag', 'var') || isempty(qr_flag))
    qr_flag = 0;
end
if(qr_flag)
    [Qc, ~] = qr(Bc, 0);
    [Qr, ~] = qr(Br', 0); % Qr=Qr'; 
else
    [Qc, ~, ~] = lansvd(Bc,r,'L','OPTIONS');
    [~, ~,Qr] = lansvd(Br,r,'L','OPTIONS');
end


B=(Ac'*Qr\Bc'*Qc);
[U S V] = svd(B');
U = Qc*U; V=Qr*V;
Xc = U*S*V';

B=(Ar*Qc\Br*Qr);
[U S V] = svd(B);
U = Qc*U; V=Qr*V;
Xr = U*S*V';
cost_Xc=norm(Ar*Xc-Br,'fro')+norm(Xc*Ac-Bc,'fro');
cost_Xr=norm(Ar*Xr-Br,'fro')+norm(Xr*Ac-Bc,'fro');

% Take solution minimizing loss function
if(cost_Xc<cost_Xr)
    X=Xc;
else
    X=Xr;
end

end

