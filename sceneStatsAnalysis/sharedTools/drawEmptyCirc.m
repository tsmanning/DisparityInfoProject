function drawEmptyCirc(x,y,r,linecolor,linewidth)

% Draw a dot of a given size, location, and color

numVerts = 50;

plot(r*cos(linspace(0,2*pi,numVerts))+x,...
     r*sin(linspace(0,2*pi,numVerts))+y,...
     'color',linecolor,'linewidth',linewidth);

end