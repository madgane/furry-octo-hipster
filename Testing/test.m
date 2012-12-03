

clc;

Pmax = 1;
H = complex(randn(1,4,4),randn(1,4,4));

cvx_begin

variable t
variables d1 d2
variable X(4,1) complex

minimize t

norm(X * H(:,:,1)) <= t


norm(X) 0= Pmax

cvx_end