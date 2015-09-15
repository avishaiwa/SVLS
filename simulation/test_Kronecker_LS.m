% Test kronecker least-squares problem 
m = 12; p = 6; A1 = randn(m, p); 
n = 6; q = 4; A2 = randn(n, q); 
x = randn(p*q, 1); 
sigma = 0.12; 

b = kron(A1, A2) * x + randn(m*n,1)*sigma; 

x_hat_kron  = Kronecker_LS(A1, A2, b)
x_hat  = Kronecker_LS(A1, A2, b, 'brute_force')

figure; plot(x_hat, x_hat_kron, '.')


figure; plot(x, x_hat, 'o'); 
hold on; plot(x, x_hat_kron, 'r*'); 
plot(x, x, 'k'); 
xlabel('true x'); ylabel('x-hat'); 
legend({'brute-force', 'kronecker'})
