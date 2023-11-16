% Définition des paramètres
a = 0;
b = 2;
beta = 1;
M = 49;
k = 0.0001;
T = 1;
N = T/k;
alpha = 1;
% Calcul du pas h
h = (b - a) / (M + 1);

% Initialisation de x et U(0)
x = linspace(a + h, b - h, M)';

U0 = zeros(M, 1);  % Initialisation de U0

% Définition de la fonction u0(x)
for i = 1:M
    if x(i) > a && x(i) <= (b + a) / 2
        U0(i) = 2 * x(i);
    else 
        U0(i) = 2 * (a + b - x(i));
    end
end

% Préallocation de U pour stocker l'historique des solutions
U = zeros(M, N);
U(:,1) = U0; % La première colonne est U0

p = 1-2*alpha*k/(h^2);
q = alpha*k/(h^2);
r = q;

% Création d'une matrice carrée MxM remplie de zéros
A = zeros(M, M);

% Remplissage de la diagonale principale avec p
for i = 1:M
    A(i, i) = p;
end

% Remplissage de la diagonale supérieure avec q
for i = 1:M-1
    A(i, i+1) = q;
end

% Remplissage de la diagonale inférieure avec r
for i = 2:M
    A(i, i-1) = r;
end


% Initialisation de la figure pour l'affichage dynamique
figure;
% Tracer la solution numérique initiale
hPlot = plot(x, U(:, 1), 'LineWidth', 2);
hold on; % Maintenir la figure pour les tracés additionnels
ylim([-0.2 2]);

for n = 2:N

    % La colonne précédente représente UOLD
    UOLD = U(:, n-1);

    % Calcul de UNEW en utilisant la matrice A et UOLD
    UNEW = A * UOLD;

    % Mise à jour de la matrice U avec les nouvelles valeurs
    U(:, n) = UNEW;

    % Mise à jour du tracé de la solution numérique
    set(hPlot, 'YData', U(:, n));

    drawnow;
    pause(0.0002);
end

hold off; % Fin du maintien de la figure