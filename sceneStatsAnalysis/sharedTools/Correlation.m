function [r,varargout] = Correlation(x,y)
% function [r[,p,t]] = Correlation(x,y)
% Calculates the correlation between x and y and its significance.
% Returns r, the product-moment correlation coefficient, Armitage & Berry p. 163.
% Optionally returns p (significance) and t.
%
% Example: 
% x = rand(1,20); y = x + rand(1,20);
% r = Correlation(x,y): returns r, the product-moment correlation coefficient.
% [r,p] = Correlation(x,y): returns r, the product-moment correlation coefficient, and also p, the probability of getting this r under the null hypothesis.
% [r,p,t] = Correlation(x,y): returns r and and also the t statistic.
%
% Jenny Read 2/4/2003

xbar = mean(x);
ybar = mean(y);

if isnan(xbar) | isnan(ybar)
    r = NaN;
    if nargout>1
        varargout{1}=NaN;
    end
    if nargout>2
        varargout{2}=NaN;
    end
    return
end

n=length(x);
r = sum( (x-xbar).*(y-ybar)) ./ sqrt( sum((x-xbar).^2) * sum((y-ybar).^2) );
if nargout>1
    t = r*sqrt((n-2)/(1-r^2));
    % work out probability of getting this t by chance
    prob = tcdf(t,n-2);
    if prob>0.5
        p = 2*(1-prob);
    else
        p = 2*prob;
    end
    if nargout>1
        varargout{1}=p;
    end
    if nargout>2
        varargout{2}=t;
    end
end
