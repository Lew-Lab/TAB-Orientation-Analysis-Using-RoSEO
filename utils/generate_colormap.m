function [map] = generate_colormap(middlecolor)

map = NaN(510,3);
for i = 1:255
    map(i,:) = middlecolor/255*i;
end
for i = 256:510
    map(i,:) = ([1,1,1]-middlecolor)/255*(i-255)+middlecolor;
end

end

