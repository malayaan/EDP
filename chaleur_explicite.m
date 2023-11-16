% Définition des paramètres
a = 0;
b = 1;
M = 49;
T = 0.001;
k_values = [0.0001, 0.0002, 0.00021, 0.00022];
alpha = 1;
h = (b - a) / (M + 1);
nbTermes = 5;
tps=0.02;
% Initialisation de x
x = linspace(a + h, b - h, M)';

% Définition de la fonction u0(x)
U0 = zeros(M, 1);
for i = 1:M
    mu_m = rand; % Génération d'une valeur aléatoire pour mu_m
    if x(i) > a && x(i) <= (b + a) / 2
        U0(i) = 2 * x(i);
    else 
        U0(i) = 2 * (a + b - x(i));
    end
    %U0(i) = - mu_m * (x(i) - a) * (x(i) - b);
    %attention signe moins a supprimer c'est juste pour l'affichage
end

for k = k_values
    N = round(T/k); % Arrondir N à l'entier le plus proche
    p = 1-2*alpha*k/(h^2);
    q = alpha*k/(h^2);
    r = q;
    fprintf('Pour k = %f, alpha*k/h^2 = %f\n', k, alpha*k/h^2);

    A = diag(p*ones(1,M)) + diag(q*ones(1,M-1),1) + diag(r*ones(1,M-1),-1);
    
    U = zeros(M, N);
    U(:,1) = U0;

    figure;
    hPlotNum = plot(x, U(:, 1), 'LineWidth', 2, 'Color', 'blue');
    hold on;
    ylim([-0.2 2]);
    title(sprintf('Solution numérique pour k = %f', k));

    if k == k_values(1)
        % Calcul initial de la solution exacte
        u_exact = calculeSolutionExacte(x, 0, alpha, nbTermes);
        hPlotExact = plot(x, u_exact, 'LineWidth', 2, 'Color', 'red');
        legend('Solution numérique', 'Solution exacte');
    else
        legend('Solution numérique');
    end

    for n = 2:N
        UOLD = U(:, n-1);
        UNEW = A * UOLD;
        U(:, n) = UNEW;

        % Mise à jour de la solution numérique
        set(hPlotNum, 'YData', U(:, n));

        if k == k_values(1)
            % Calcul et mise à jour de la solution exacte
            t = (n-1) * k;
            u_exact = calculeSolutionExacte(x, t, alpha, nbTermes);
            set(hPlotExact, 'YData', u_exact);
        end

        drawnow;
        pause(tps);
    end
    pause(1)
    hold off;
    if k == k_values(1)
        
    end
end

function u_exact = calculeSolutionExacte(x, t, alpha, nbTermes)
    u_exact = zeros(size(x));
    for j = 0:(nbTermes - 1)
        terme = ((-1)^j / ((2*j + 1)^2 * pi^2)) * sin((2*j + 1) * pi * x) * exp(-(2*j + 1)^2 * pi^2 * alpha * t);
        u_exact = u_exact + terme;
    end
    u_exact = u_exact * 8;
end
