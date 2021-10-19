function k = image_to_index(x,y,d)

if x~=d
    k = d*(y-1)+x;
else
    k = d*y;
end