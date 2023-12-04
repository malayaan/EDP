% Définition des paramètres
a = 0;
b = 2;
beta = 1/2;
M = 200;
k = 0.005;
T = 2;
N = T/k;
theta=1;

% Calcul du pas h
h = (b - a) / (M + 1);

% Initialisation de x et U(0)
x = linspace(a + h, b-h, M+1)';
U0 = exp(-5 * (5 * x - 1) .^ 2);

% Préallocation de U pour stocker l'historique des solutions
U = zeros(M+1, N);
U(:,1) = U0; % La première colonne est U0

p = 1/k-2*theta;
q = theta;
r = theta;

% Création d'une matrice carrée MxM remplie de zéros
A = zeros(M+1, M+1);

% Remplissage de la diagonale principale avec p
for i = 1:M+1
    if i ~= M+1
        A(i, i) = 1/k-2*theta;
    else
        A(i, i) = 1;
    end
end

% Remplissage de la diagonale supérieure avec q
for i = 1:M
    if i ~= M+1
        A(i, i+1) = theta;
    else
        A(i, i+1) = 0;
    end
end

% Remplissage de la diagonale inférieure avec r
for i = 2:M+1
    if i ~= M+1
        A(i, i-1) = theta;
    else
        A(i, i-1) = 0;
    end
end

% Création d'une matrice carrée MxM remplie de zéros
B = zeros(M+1, M+1);

% Remplissage de la diagonale principale avec p
for i = 1:M+1
    if i ~= M+1
        B(i, i) = 1/k-2*beta/h^2*(1-theta);
    else
        B(i, i) = 1-beta*k/h;
    end
end

% Remplissage de la diagonale supérieure avec q
for i = 1:M
    if i ~= M+1
        B(i, i+1) = beta/h^2*(1-theta);
    else
        B(i, i+1) =0;
    end
end

% Remplissage de la diagonale inférieure avec r
for i = 2:M+1
    if i ~= M+1
        B(i, i-1) = beta/h^2*(1-theta);
    else
        B(i, i-1) = beta * k / h;
    end
end

C=B/A;

% Initialisation de la figure pour l'affichage dynamique
figure;
% Tracer la solution numérique initiale
hPlot = plot(x, U(:, 1), 'LineWidth', 2);
hold on; % Maintenir la figure pour les tracés additionnels
% Tracer la solution exacte initiale (sera mise à jour dans la boucle)
hExactPlot = plot(x, U0, 'r--', 'LineWidth', 2); % en rouge pointillé pour la solution exacte
ylim([-0.2 1]);

for n = 2:N
    % La colonne précédente représente UOLD
    UOLD = U(:, n-1);

    % Calcul de UNEW en utilisant la matrice A et UOLD
    UNEW = C * [UOLD];

    % Mise à jour de la matrice U avec les nouvelles valeurs
    U(:, n) = UNEW(1:M+1);

    % Mise à jour du tracé de la solution numérique
    set(hPlot, 'YData', U(:, n));

    % Calcul de la solution exacte à l'instant t
    xExact = x - beta*k*(n-1);
    UExact = exp(-5 * (5 * xExact - 1) .^ 2);

    % Mise à jour du tracé de la solution exacte
    set(hExactPlot, 'XData', x, 'YData', UExact);

    drawnow;
    pause(0.002);
end

legend([hPlot, hExactPlot], 'Solution numérique', 'Solution exacte'); % Ajouter une légende
hold off; % Fin du maintien de la figure