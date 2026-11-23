function Hq = cubicsplineinterp(X, H, Xq)
% Interpolates through cubic spline.
% INPUTS:
%   X  - Independent variable data points (Temp or Pressure)
%   H  - Independent corresponding enthalpies
%   Xq - Evaluation point(s)
% OUTPUT:
%   Hq - Interpolated enthalpy value(s)


% Convert input data into column vectors
X = X(:);
H = H(:);

% Validate data; at least two points, increasing, matching sizes so data points correspond
if length(X) ~= length(H)
    error('X and H must have the same number of elements.');
end

if length(X) < 2
    error('At least two data points are required.');
end

if any(diff(X) <= 0)
    error('X values must be strictly increasing.');
end

% Check evaluation points in case they are not normal numerical values
if ~isnumeric(Xq) || ~isreal(Xq) || isempty(Xq) || any(~isfinite(Xq(:)))
    error('Evaluation points must be finite real values.');
end

% Prevent extrapolation outside of the provided values
if any(Xq(:) < X(1)) || any(Xq(:) > X(end))
    error('Evaluation point is outside the data range.');
end

% Number of inputs
n = length(X);

% Interval widths (step size should be consistent for our use but keeping code usable if it isn't)
dx = diff(X);

% Empty matrices
A = zeros(n,n);
B = zeros(n,1);

% Matrix boundaries
A(1,1) = 1;
A(n,n) = 1; 

% Remaining matrix values
for i = 2:n-1
    A(i,i-1) = dx(i-1);
    A(i,i) = 2*(dx(i-1)+dx(i));
    A(i,i+1) = dx(i);
    B(i) = 6*((H(i+1)-H(i))/dx(i) - (H(i)-H(i-1))/dx(i-1));
end

% Second derivatives
M = A\B;

% Interval with desired point
i = find(Xq >= X(1:end-1) & Xq <= X(2:end));

% Cubic Coefficients
a = (H(i+1)-H(i))/dx(i) - dx(i)*(2*M(i)+M(i+1))/6;
b = M(i)/2;
c = (M(i+1)-M(i))/(6*dx(i));

% Final Cubic Equation
tp = Xq-X(i);

Hq = H(i) + a*tp + b*tp^2 + c*tp^3;

end

% Remove comment on the following for easy test case - X is temp in...
% celsius, H is enthalpy at 0.05 MPa - tabulated Hq is 2731.5
% X = [100, 110, 120, 130, 140, 150];
% H = [2682.4, 2702.1, 2721.7, 2741.2, 2760.7, 2780.2];
% Xq = 125;

% Hq = cubicsplineinterp(X, H, Xq);
% fprintf('Enthalpy: %f\n', Hq)
