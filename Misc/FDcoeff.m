%% Finite Difference Formulas (Fornberg 1988)

M = 3;      % Order of the highest derivative
N = 4;      % Number of grid points - 1

x0 = 0;     % Extrapolation point

% Grid points
alpha = (0:N) + 4.5;

%% Calculate coefficients

% Initialize array
delta = zeros(M+1, N+1, N+1);
delta(1,1,1) = 1.0;
c1 = 1.0;

% Algorithm
for n = 1:N
    c2 = 1.0;

    for nu = 0:n-1
        c3 = alpha(n+1) - alpha(nu+1);
        c2 = c2 * c3;

        % m = 0
        delta(1,n+1,nu+1) = (alpha(n+1)-x0) * delta(1,n,nu+1) / c3;

        for m = 1:min(n,M)
            delta(m+1,n+1,nu+1) = ( (alpha(n+1)-x0) * delta(m+1,n,nu+1) ...
                - m*delta(m,n,nu+1) ) / c3;
        end
    end

    % m = 0
    delta(1,n+1,n+1) = c1/c2 * ( - (alpha(n)-x0) * delta(1,n,n) );

    for m = 1:min(n,M)
        delta(m+1,n+1,n+1) = c1/c2 * ( m*delta(m,n,n) ...
            - (alpha(n)-x0) * delta(m+1,n,n) );
    end
    c1 = c2;
end

%% Display coefficients

m = 0;      % Derivative order
n = 2;      % Grid points - 1

fprintf('Deriv. order: %d\n', m);
fprintf('Extrapo. x0 = %f\n', x0);
fprintf('Grid point coords:\n');
disp(alpha(1:n+1));
fprintf('Coefficient delta:\n');
disp(squeeze(delta(m+1, n+1, 1:n+1))');
fprintf('var = %f\n', sum(alpha(1:n+1) .* squeeze(delta(m+1, n+1, 1:n+1))'));

