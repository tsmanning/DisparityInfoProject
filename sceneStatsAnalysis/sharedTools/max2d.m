function [mx,varargout] = max2d(arr)
% Finds the maximum of a 2d array, and (optionally) outputs the indices of
% the maximum.
% Eg, suppose z is an array of size (ny,nx) (so you can do imagesc(x,y,z)).
% mx = max2d(z) outputs max(max(z)), the overall maximum of z.
% [mx,jx,jy] = max2d(z) also gives you jx,jy such that z(jy,jx) = mx, ie
% the position of the peak is at x(jx) and y(jy).


[ny,nx] = size(arr);
arr = reshape(arr,[1 nx*ny]);
[mx,jmx] = max(arr);
if nargout>1
    jx = ceil(jmx/ny);
    varargout{1}=jx;
    if nargout>2
        if jmx/ny==fix(jmx/ny)
            jy = ny;
        else
            jy = mod(jmx,ny);
        end
        varargout{2}=jy;
    end
end

