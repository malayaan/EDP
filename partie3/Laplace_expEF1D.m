% Définition des bornes de l'intervalle
a = 0; % Borne inférieure de l'intervalle
b = 1; % Borne supérieure de l'intervalle

n = 5; % Nombre de points intérieurs
h = (b - a) / (n + 1); % Pas de découpage

%initialisation du temps d'éxécution et du nombre de points temporels
T = 2;
m = 20;

% Initialisation des matrices P1 et P2
P2 = zeros(n, n);

% Définitions des fonctions psi et dpsi
psi = @(i, x) (x - (a + (i-1)*h)) / h .* ((a + (i-1)*h) <= x & x <= (a + i*h)) + ...
               ((a + (i+1)*h) - x) / h .* ((a + i*h) <= x & x <= (a + (i+1)*h));
dpsi = @(i, x) (1 / h) * ((a + (i-1)*h) <= x & x <= (a + i*h)) - ...
                (1 / h) * ((a + i*h) <= x & x <= (a + (i+1)*h));

% Calcul de P2
for i = 1:n
    % Éléments diagonaux
    f = @(x) psi(i, x).^2 / h;
    P2(i, i) = integral(f, a + (i-1)*h, a + (i+1)*h);
    
    % Éléments hors diagonale
    if i < n
        g = @(x) psi(i, x).*psi(i+1, x) / h;
        value = integral(g, a + i*h, a + (i+1)*h);
        P2(i, i+1) = value;
        P2(i+1, i) = value; % La matrice est symétrique
    end
end

% P1 peut maintenant être calculé à partir de P2 en ajoutant les termes dérivés
P1 = P2; % Commencer avec les valeurs de P2
for i = 1:n
    % Ajouter les termes dérivés aux éléments diagonaux
    f_prime = @(x) dpsi(i, x).^2;
    P1(i, i) = P1(i, i) + integral(f_prime, a + (i-1)*h, a + (i+1)*h);
    
    % Ajouter les termes dérivés aux éléments hors diagonale
    if i < n
        g_prime = @(x) dpsi(i, x).*dpsi(i+1, x);
        value = integral(g_prime, a + i*h, a + (i+1)*h);
        P1(i, i+1) = P1(i, i+1) + value;
        P1(i+1, i) = P1(i+1, i) + value; % La matrice est symétrique
    end
end

% Affichage des matrices
disp(P2);
disp(P1);

% Points de maillage
x = linspace(a, b, n+2); % n+2 points incluant les bornes

% Fonction initiale
u_init = @(x) exp(-5 * (5 * x - 1) .^ 2);

% Évaluation de la fonction initiale aux points de maillage intérieurs
A0 = u_init(x(2:end-1)); % Exclure les points aux bornes

% Affichage du vecteur A0
%disp(A0);

% Définition de la fonction f(t,x)
f = @(t, x) 0;

% Définition de la fonction B(t)
B = @(t) arrayfun(@(j) integral(@(x) f(t, x) .* psi(j, x), a, b), 1:n);

% Pas de temps pour la discrétisation temporelle
dt = T / m;

% Préallocation de A pour stocker l'historique des solutions
A = zeros(n, m);
A(:, 1) = A0; % La première colonne est A0

% Initialisation de la figure pour l'affichage dynamique
figure;
% Tracer la solution numérique initiale
hPlot = plot(x(2:end-1), A(:, 1), 'LineWidth', 2);
hold on; % Maintenir la figure pour les tracés additionnels
ylim([-0.5, 1]); % Limite pour l'axe des ordonnées
xlim([a, b]); % Limite pour l'axe des abscisses


for M = 2:m
    % La colonne précédente représente AOLD
    AOLD = A(:, M-1);
    
    % Calcul de ANEW 
    ANEW = P2 \ (B(dt*(M))' + P1*AOLD); % Utilisation de l'opérateur de division matricielle à gauche pour résoudre P1*ANEW = ...
    
    % Mise à jour de la matrice A avec les nouvelles valeurs
    A(:, M) = ANEW;

    % Mise à jour du tracé de la solution numérique
    set(hPlot, 'YData', A(:, M));

    drawnow;
    pause(0.2); % Pause pour permettre l'affichage du tracé
end

hold off; % Fin du maintien de la figure
