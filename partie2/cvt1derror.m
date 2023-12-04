% Définition des paramètres
a = 0;
b = 2;
beta = 1;
M_values = [49, 99, 199, 399, 799]; % Différentes valeurs de M
T = 1.0;

% Configuration du graphique
figure;
hold on;
colors = ['b', 'g', 'r', 'c', 'm']; % Couleurs pour différentes valeurs de M
legendInfo = cell(1, length(M_values) * 2); % Pour la légende

for idx = 1:length(M_values)
    M = M_values(idx);
    h = (b - a) / (M + 1);
    k = 0.5 * h; % Définir k pour que k/h = 0.5
    N = round(T/k);
    x = linspace(a + h, b - h, M)';
    U0 = exp(-5 * (5 * x - 1) .^ 2);

    U = zeros(M, N);
    U(:,1) = U0;

    for n = 1:(N-1)
        for m = 2:M-1
            U(m, n+1) = (1 - k*beta/h) * U(m, n) + (k * beta / h) * U(m-1, n);
        end

        % Solution exacte à l'instant t
        xExact = x - beta*k*n;
        UExact = exp(-5 * (5 * xExact - 1) .^ 2);
    end

    % Tracé de la solution numérique
    plot(x, U(:, end), 'Color', colors(idx), 'LineWidth', 2);
    legendInfo{idx * 2 - 1} = ['Numerical, M = ', num2str(M)];

    % Tracé de la solution exacte
    plot(x, UExact, [colors(idx) '--'], 'LineWidth', 2);
    legendInfo{idx * 2} = ['Exact, M = ', num2str(M)];

    % Calcul de l'erreur L2
    L2_error = norm(U(:, end) - UExact, 2);
    disp(['L2 Error for M = ' num2str(M) ': ' num2str(L2_error)]);

    % Trouver le maximum de la solution numérique et sa position
    [maxU, pos] = max(U(:, end));
    disp(['Maximum value for M = ' num2str(M) ': ' num2str(maxU) ' at position ' num2str(x(pos))]);
end

title('Comparaison des solutions numériques et exactes pour différentes valeurs de M');
xlabel('x');
ylabel('U(x,T)');
legend(legendInfo);
hold off;
