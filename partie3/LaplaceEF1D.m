a = 0;
b = 1;
n = 100; % Nombre de points intérieurs
h = (b - a) / (n + 1); % Le pas h

X = linspace(a, b, n+2); % Vecteur contenant n+2 points, y compris a et b

% Initialisation des fonctions psi
psi = cell(1, n); % Cell array pour stocker les fonctions pour les n points intérieurs

% Création des fonctions de base psi pour les points intérieurs
for i = 2:n+1
    psi{i} = @(X_query) ((X_query >= X(i-1) & X_query < X(i)) .* (X_query - X(i-1)) / h) +  ((X_query >= X(i) & X_query < X(i+1)) .* (X(i+1) - X_query) / h);
end

%definition de f
f = @(x) pi^2 * sin(pi * x);

% Création d'une matrice carrée n*n remplie de zéros
M = zeros(n, n);
p=2/h;
q=-1/h;

% Remplissage de la diagonale principale avec p
for i = 1:n
    M(i, i) = p;
end

% Remplissage de la diagonale supérieure avec q
for i = 1:n-1
    M(i, i+1) = q;
end

% Remplissage de la diagonale inférieure avec q
for i = 2:n
    M(i, i-1) = q;
end

% Création d'une matrice colonne n*1 remplie de zéros
B = zeros(n,1);

% Remplissage de B en utilisant la méthode des trapèzes pour les n points intérieurs
for j = 2:n+1
    B(j-1) = h/2 * (f(X(j-1)) * psi{j}(X(j-1)) + 2 * f(X(j)) * psi{j}(X(j)) + f(X(j+1)) * psi{j}(X(j+1)));
end


A_solution = inv(M) * B;

% Calcul de la fonction u sur la grille de points X
u = zeros(size(X));

% Ajout des contributions de chaque fonction de base psi_i pondérée par son coefficient dans A_solution
for i = 2:n+1
    u = u + A_solution(i-1) * psi{i}(X);
end

% Tracé de la fonction u
figure;
plot(X, u, 'LineWidth', 2);
title('Fonction u approximée par la méthode des différences finies');
xlabel('x');
ylabel('u(x)');

% Calcul de sin(pi*x) sur la grille de points X
Y = sin(pi * X);

% Tracé de la fonction sin(pi*x)
hold on; % Garde le tracé précédent sur la figure
plot(X, Y, 'r--', 'LineWidth', 2); % Tracer en rouge pointillé
legend('u approximée', 'sin(pi*x)'); % Ajouter une légende
hold off; % Fin du mode hold on