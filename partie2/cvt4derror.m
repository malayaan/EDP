% Définition des paramètres
a = 0;
b = 2;
beta = 1;
T = 1;
M_values = [49, 99, 199, 399, 799]; % Différentes valeurs de M

% Initialisation de la figure
figure;
hold on;

% Boucle sur différentes valeurs de M
for M = M_values
    % Calcul de h et ajustement de k
    h = (b - a) / (M + 1);
    k = 0.5 * h;
    N = round(T / k);

    % Initialisation de x et U(0)
    x = linspace(a + h, b - h, M)';
    U0 = exp(-5 * (5 * x - 1) .^ 2);

    % Initialisation et préallocation de U
    U = zeros(M, N);
    U(:, 1) = U0;

    % Définition des coefficients p et q
    p = -k * beta / h;
    q = -p;

    % Construction de la matrice A
    A = diag(p * ones(1, M-1), 1) + diag(q * ones(1, M-1), -1);

    % Calcul de la deuxième colonne U1
    U1 = zeros(M, 1);
    for m = 2:M-1
        U1(m) = U0(m) + 0.5 * k * (beta / h * (U0(m+1) - U0(m-1)));
    end
    U(:, 2) = U1;

    % Boucle sur le temps
    for n = 3:N
        % Calcul de UNEW
        UNEW = A * U(:, n-1) + U(:, n-2);

        % Mise à jour de la matrice U
        U(:, n) = UNEW;
    end

    % Tracé de la solution numérique pour M actuel
    plot(x, U(:, end), 'DisplayName', ['M = ' num2str(M)]);

    % Calcul de la solution exacte
    xExact = x - beta * T;
    UExact = exp(-5 * (5 * xExact - 1) .^ 2);

    % Calcul de l'erreur L2 et du maximum
    L2_error = norm(U(:, end) - UExact, 2);
    [maxU, maxPos] = max(U(:, end));
    disp(['M = ' num2str(M) ', Max U: ' num2str(maxU) ', Position: ' num2str(x(maxPos)) ', L2 Error: ' num2str(L2_error)]);
end

% Tracé de la solution exacte
plot(x, UExact, 'r--', 'DisplayName', 'Solution exacte');
title('Comparaison des solutions numériques et exactes pour différentes valeurs de M');
xlabel('x');
ylabel('U(x, T)');
legend;
hold off;
