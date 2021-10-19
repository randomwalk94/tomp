function [x,y] = index_to_image(k,d)
k = k+d;
y = floor(k/d);
x = k-y*d;
if x==0
    x = d;
    y = y-1;
end