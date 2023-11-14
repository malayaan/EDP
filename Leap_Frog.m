% Définition des paramètres
a = 0;
b = 2;
beta = 1;
M = 99;
k = 0.02;
T = 1;
N = T/k;

% Calcul du pas h
h = (b - a) / (M + 1);

% Initialisation de x et U(0)
x = linspace(a + h, b - h, M)';
U0 = exp(-5 * (5 * x - 1) .^ 2);
U1 = zeros(M, 1);
for m = 2:M-1
    U1(m) = (1 - k*beta/h) * U0(m) + (k * beta / h) * U0(m-1);
end

% Préallocation de U pour stocker l'historique des solutions
U = zeros(M, N);
U(:,1) = U0; % La première colonne est U0
U(:,2) = U1; % La deuxieme colonne est U1


p = -k*beta/h;
q = - p;

% Création d'une matrice carrée MxM remplie de zéros
A = zeros(M, M);

% Remplissage de la diagonale supérieure avec p
for i = 1:M-1
    A(i, i+1) = p;
end

% Remplissage de la diagonale inférieure avec q
for i = 2:M
    A(i, i-1) = q;
end

%initialisation de la condition au limite
Uborder = zeros(N, 1);
Uborder(1,1) = exp(-5 * (5 * (b) - 1) .^ 2);
Uborder(2,1) = -beta*k/h*(Uborder(1, 1) - U(M, 1)) + k*Uborder(1, 1);

% Initialisation de la figure pour l'affichage dynamique
figure;
% Tracer la solution numérique initiale
hPlot = plot(x, U(:, 1), 'LineWidth', 2);
hold on; % Maintenir la figure pour les tracés additionnels
% Tracer la solution exacte initiale (sera mise à jour dans la boucle)
hExactPlot = plot(x, U0, 'r--', 'LineWidth', 2); % en rouge pointillé pour la solution exacte
ylim([-0.2 1]);

for n = 3:N
    Uborder(n, 1) = -beta*k/h*(Uborder(n-1, 1) - U(M, n-1)) + k*Uborder(n-1, 1);

    % La colonne précédente représente UOLD
    UOLD = U(:, n-1);
    UOLDOLD = U(:, n-2);
    % Calcul de UNEW en utilisant la matrice A et UOLD
    UNEW = A * UOLD+UOLDOLD;
    UNEW(M) = UNEW(M) + p*Uborder(n-1, 1);

    % Mise à jour de la matrice U avec les nouvelles valeurs
    U(:, n) = UNEW;

    % Mise à jour du tracé de la solution numérique
    set(hPlot, 'YData', U(:, n));

    % Calcul de la solution exacte à l'instant t
    xExact = x - beta*k*(n-1);
    UExact = exp(-5 * (5 * xExact - 1) .^ 2);

    % Mise à jour du tracé de la solution exacte
    set(hExactPlot, 'XData', x, 'YData', UExact);

    drawnow;
    pause(0.02);
end

legend([hPlot, hExactPlot], 'Solution numérique', 'Solution exacte'); % Ajouter une légende
hold off; % Fin du maintien de la figure