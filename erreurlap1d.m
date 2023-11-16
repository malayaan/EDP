% Cas à traiter
cas = {'a', 'b', 'c'};

for i = 1:length(cas)
    fprintf('Cas %s:\n', cas{i});
    fprintf('%-8s %-8s %-8s %-8s\n', 'M', 'eM', 'eM/h', 'eM/h^2');

    for K = 2:11
        M = 2 ^ K;

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
        
        % Calcul de l'erreur
        eM = max(abs(U - exact));
        
        % Affichage des résultats
        fprintf('%-8d %-8.4f %-8.4f %-8.4f\n', M, eM, eM/h, eM/h^2);
    end
    fprintf('\n');
end
