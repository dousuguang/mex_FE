clc
clear

if 1
    mex element.c;
end

tic
for i=1:10000
    [k,m]=element(1,1,1,1,1,[2,555],1,1);
end
toc

tic
for i=1:10000
    [k,m]=element_m(1,1,1,1,1,[2,555],1,1);
end
toc