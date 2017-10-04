function x = tridiag( a, b, c, f )

n = length(f);
v = zeros(n,1);   
x = v;
w = b(1);
x(1) = f(1)/w;
for i=2:n
    v(i-1) = c(i-1)/w;
    w = b(i) - a(i)*v(i-1);
    x(i) = ( f(i) - a(i)*x(i-1) )/w;
end
for j=n-1:-1:1
   x(j) = x(j) - v(j)*x(j+1);
end


%% Credits to: Mark Holmes