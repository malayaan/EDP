% Paramètres du problème
a = 0;
b = 1;
h_values = [1/10, 1/20, 1/40, 1/80]; % Différents pas h pour l'analyse de convergence

% Initialisation des tableaux pour stocker les erreurs et les pas h
errors_L2 = zeros(length(h_values), 1);
h_log = log(h_values);

%definition de f
f = @(x) pi^2 * sin(pi * x);

% Boucle sur les différents pas h
for idx = 1:length(h_values)
    h = h_values(idx);
    n = (b - a) / h - 1; % Nombre de points intérieurs
    X = linspace(a, b, n+2); % Vecteur contenant n+2 points, y compris a et b
    
    % Initialisation des fonctions psi
    psi = cell(1, n); % Cell array pour stocker les fonctions pour les n points intérieurs

    % Création des fonctions de base psi pour les points intérieurs
    for i = 2:n+1
        psi{i} = @(X_query) ((X_query >= X(i-1) & X_query < X(i)) .* (X_query - X(i-1)) / h) +  ((X_query >= X(i) & X_query < X(i+1)) .* (X(i+1) - X_query) / h);
    end

    % Matrice M et vecteur B
    M = (2/h) * diag(ones(n, 1)) - (1/h) * diag(ones(n-1, 1), 1) - (1/h) * diag(ones(n-1, 1), -1);
    % Création d'une matrice colonne n*1 remplie de zéros
    B = zeros(n,1);
    % Remplissage de B en utilisant la méthode des trapèzes pour les n points intérieurs
    for j = 2:n+1
        B(j-1) = h/2 * (f(X(j-1)) * psi{j}(X(j-1)) + 2 * f(X(j)) * psi{j}(X(j)) + f(X(j+1)) * psi{j}(X(j+1)));
    end
    
    
    A_solution = inv(M) * B;
    
    x = linspace(a + h, b - h, n); % Points intérieurs seulement
    U_exact = sin(pi * x);

    % Calcul de l'erreur e = U_exact - U_approx
    e = U_exact - A_solution';
    %disp(e)

    % Matrice de masse W (simplifiée pour des fonctions linéaires par morceaux)
    W = (h/6) * (diag(ones(n-1, 1), -1) + 4*diag(ones(n, 1), 0) + diag(ones(n-1, 1), 1));
    
    % Calcul de la norme L2 de l'erreur en utilisant la matrice de masse
    errors_L2(idx) = sqrt(e * W * e');

    % Affichage de la norme L2 de l'erreur pour le pas h courant
    fprintf('Norme L2 de lerreur pour h = %f : %f\n', h, errors_L2(idx));
end

% Tracer la courbe de convergence en échelle log-log
error_log = log(errors_L2);
plot(h_log, error_log, '-o');
xlabel('log(h)');
ylabel('log(Error L2)');

% Ajustement linéaire pour estimer la pente
fit_params = polyfit(h_log, error_log, 1);
slope = fit_params(1);

% Affichage de l'ordre de convergence estimé
fprintf('Lordre estimé de convergence en norme L2 est : %f\n', -slope);
