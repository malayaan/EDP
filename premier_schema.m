%on définit les valeurs initiales
% Définition des paramètres
a = 0;
b = 2;
beta = 1;
M = 99;
k = 0.01;
T = 1;
N = T/k;

% Calcul du pas h
h = (b - a) / (M + 1);

% Initialisation de x et U(0)
x = linspace(a + h, b - h, M)';
U0 = exp(-5 * (5 * x - 1) .^ 2);

% Préallocation de U pour stocker l'historique des solutions
U = zeros(M, N);
U(:,1) = U0; % La première colonne est U0

% Configuration du graphique pour l'affichage dynamique
figure;
hPlot = plot(x, U0, 'LineWidth', 2);
hold on; % Maintient la figure actuelle pour le tracé de la solution exacte
hExactPlot = plot(x, U0, 'r--', 'LineWidth', 2); % Solution exacte en rouge pointillé
axis([a b -0.1 1.1]); % Réglez les axes selon vos données et préférences
title('Évolution de U avec le temps');
xlabel('x');
ylabel('U(x,t)');
drawnow;

% Boucle pour le calcul itératif de U(n+1) avec affichage dynamique
for n = 1:N-1
    % Mise à jour de la solution numérique U
    for m = 2:M-1 % Commencer à m = 2 et finir à M-1 pour ne pas dépasser les limites de la matrice
        U(m, n+1) = (1 - k*beta/h) * U(m, n) + (k * beta / h) * U(m-1, n);
    end
    ù
    % Mise à jour du graphique pour la solution numérique
    set(hPlot, 'YData', U(:, n+1));
    
    % Calcul de la solution exacte à l'instant t
    xExact = x - beta*k*n;
    UExact = exp(-5 * (5 * xExact - 1) .^ 2);
    
    % Mise à jour du graphique pour la solution exacte
    set(hExactPlot, 'XData', x, 'YData', UExact);
    
    drawnow; % Met à jour la figure avec les nouveaux tracés
    pause(0.1); % Pause pour une meilleure visualisation
end

legend([hPlot, hExactPlot], 'Numerical', 'Exact'); % Ajouter une légende après la création des tracés
hold off; % Relâche la figure pour de futurs tracésoikioo