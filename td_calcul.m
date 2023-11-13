N = 256; % Taille de la matrice

% Initialiser A avec des zéros
A = zeros(N);

% Construction de la matrice A
A = 4*eye(N) + diag(ones(N-1, 1), 1) + diag(ones(N-1, 1), -1);
A = A + diag(ones(N-N/2, 1), N/2) + diag(ones(N-N/2, 1), -N/2);

% Vérifier que A est symétrique
if issymmetric(A)
    disp('La matrice A est symétrique.');
else
    disp('La matrice A n''est pas symétrique.');
end

% Vecteur f aléatoire
f = rand(N, 1);

% Mesure de la taille mémoire de A
s = whos('A');
memorySize = s.bytes;

% (a) inv(A) * f
inverseMultTime = 0;
for i = 1:10
    tic; % Démarre le chronomètre
    y = inv(A) * f;
    inverseMultTime = inverseMultTime + toc; % Arrête le chronomètre et cumule le temps
end

% (b) A\f
solveTime = 0;
for i = 1:10
    tic; % Démarre le chronomètre
    x = A\f;
    solveTime = solveTime + toc; % Arrête le chronomètre et cumule le temps
end

% Affichage des résultats
fprintf('Taille mémoire de A : %d bytes\n', memorySize);
fprintf('Temps total pour inv(A) * f : %f seconds\n', inverseMultTime);
fprintf('Temps total pour A \\ f : %f seconds\n', solveTime);