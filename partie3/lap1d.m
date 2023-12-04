% Cas à traiter
cas = {'a', 'b', 'c'};
M = 20;

for i = 1:length(cas)
    switch cas{i}
        case 'a'
            % Cas (a)
            a = 0;
            b = pi;
            f = @(x) sin(x);
            solution_exacte = @(x) sin(x);
        case 'b'
            % Cas (b)
            a = -1;
            b = 2;
            f = @(x) (x - a) .* (x - b);
            solution_exacte = @(x) (-1/12) * x.^4 + (1/6) * x.^3 + x.^2 - 1.08333 * x - 1.83333;
        case 'c'
            % Cas (c)
            a = -sqrt(2);
            b = sqrt(2);
            f = @(x) -(x.^2 + 4*x) .* exp(x);
            solution_exacte = @(x) exp(x) .* (x.^2 - 2);
    end

    % Calculs communs
    h = (b - a) / (M + 1);
    p = 2/h^2;
    q = -1/h^2;
    r = q;
    A = diag(repmat(p, M, 1)) + diag(repmat(q, M-1, 1), 1) + diag(repmat(r, M-1, 1), -1);
    x = linspace(a + h, b - h, M)';
    
    % Calcul de U et de la solution exacte
    U = A \ f(x);
    exact = solution_exacte(x);
    
    % Affichage des courbes
    figure;
    plot(x, U, 'b');
    hold on;
    plot(x, exact, 'r--');
    hold off;
    title(['Cas ', cas{i}]);
    xlabel('x');
    ylabel('U et Solution exacte');
    legend('U approximé', 'Solution exacte');
end
