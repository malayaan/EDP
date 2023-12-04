% Paramètres
alpha = 1;
T = 0.5;
h_values = [0.1, 0.05, 0.025, 0.0125, 0.00625];
nbTermes = 1; % Nombre de termes pour la solution exacte
erreurs = zeros(size(h_values)); % Stocker les erreurs

% Boucle sur les différentes valeurs de h
for i = 1:length(h_values)
    h = h_values(i);
    k = 0.25 * h^2 / alpha; % Calcul de k
    M = round((2 - 0) / h) - 1;
    N = round(T / k);

    x = linspace(h, 2-h, M)';
    U = zeros(M, N);

    % Définition de la fonction u0(x)
    U0 = zeros(M, 1);
    for j = 1:M
        mu_m = rand; % Génération d'une valeur aléatoire pour mu_m
        if x(j) > a && x(i) <= (b + a) / 2
            U0(j) = 2 * x(j);
        else 
            U0(j) = 2 * (a + b - x(j));
        end
        %U0(i) = - mu_m * (x(i) - a) * (x(i) - b);
        %attention signe moins a supprimer c'est juste pour l'affichage
    end
    U(:,1) = U0;
    % Matrice A pour la solution numérique
    p = 1 - 2 * alpha * k / h^2;
    q = alpha * k / h^2;
    A = diag(p * ones(1, M)) + diag(q * ones(1, M-1), 1) + diag(q * ones(1, M-1), -1);

    % Calcul de la solution numérique
    for n = 2:N
        U(:, n) = A * U(:, n-1);
    end

    % Calcul de la solution exacte à T
    u_exact = calculeSolutionExacte(x, T, alpha, nbTermes);

    % Calcul de l'erreur en norme L∞
    erreur = max(abs(U(:, N) - u_exact));
    erreurs(i) = erreur;
end

% Tracer l'erreur en fonction de h
figure;
loglog(h_values, erreurs, '-o');
xlabel('h');
ylabel('Erreur L^\infty');
title('Erreur L^\infty en fonction de h');

function u_exact = calculeSolutionExacte(x, t, alpha, nbTermes)
    u_exact = zeros(size(x));
    for j = 0:(nbTermes - 1)
        terme = ((-1)^j / ((2*j + 1)^2 * pi^2)) * sin((2*j + 1) * pi * x) * exp(-(2*j + 1)^2 * pi^2 * alpha * t);
        u_exact = u_exact + terme;
    end
    u_exact = u_exact * 8;
end
