function FM = FixMap(X,A)
%FIXMAP Summary of this function goes here
%   Detailed explanation goes here

    FM = (A'*A)*X + A*X*A+A'*X*A'+X*(A'*A);
    
end
