function T = transformation(R)
% d = xj-xi; L = norm(d); c = d(1)/L; s = d(2)/L;
c = R(1); s = R(2);
T = [ 
	c   s  0   0  0  0
    -s  c  0   0  0  0
    0   0  1   0  0  0
    0   0  0   c  s  0
    0   0  0  -s  c  0
    0   0  0   0  0  1 ];
end
