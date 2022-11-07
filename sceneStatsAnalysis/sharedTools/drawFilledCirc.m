function drawFilledCirc(x,y,r,fillcolor,alphaL)

% Draw a dot of a given size, location, and color

numVerts = 50;

fill(r*cos(linspace(0,2*pi,numVerts))+x,...
     r*sin(linspace(0,2*pi,numVerts))+y,...
     fillcolor,...
     'FaceAlpha',alphaL);

end