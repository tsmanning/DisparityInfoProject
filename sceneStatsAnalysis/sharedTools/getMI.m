function [MI] = getMI(px,py,pxy)

% Calculate the mutual information between two variables given their joint
% and marginal distributions

% Assume columns of pxy are different x vals and rows are different y vals

[numY,numX] = size(pxy);

if ~iscolumn(py)
    py = py';
end

if ~isrow(px)
    px = px';
end

pxRep = repmat(px,[numY 1]);
pyRep = repmat(py,[1 numX]);

% should probably do some matrix instead of index math here
pInd = pxRep.*pyRep;

% Avoid issues with log & dividing by zero by making zeros eps
pxy(pxy==0) = eps;
pInd(pInd==0) = eps;

MI = sum(sum( pxy .* log(pxy./pInd) ));

end