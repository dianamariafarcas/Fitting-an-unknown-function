%% Clear and load the data
clear
load('proj_fit_06.mat');
clc

X1 = id.X{1};
X2 = id.X{2};

valX1 = val.X{1}';
valX2 = val.X{2}';

MSEs_val = zeros(25,1); % Array to store the MSEs for validation data
MSEs_id = zeros(25,1); % Array to store the MSEs for identification data
m_min = 1;  % Track the order of the best model
l_min = 3;
MSE_index = 1;

%% Convert the matrices Y into vectors

Yn_id = zeros(length(X1) * length(X1), 1);
k = 1;
for i = 1 : length(X1)
    for j = 1 : length(X2)
        Yn_id(k) = id.Y(i, j);
        k = k + 1;
    end
end

Yn_val = zeros(length(valX1) * length(valX2), 1);
k = 1;
for i = 1 : length(valX1)
    for j = 1 : length(valX2)
        Yn_val(k) = val.Y(i, j);
        k = k + 1;
    end
end

% Initial calculation
R = zeros(length(X1) * length(X2), 3); 
R = calculateR(X1, X2, 3, 1); % Calculate R matrix for identification data, for m=1
theta = R \ Yn_id; % Estimate theta
idY = R * theta; % Calculate the fitted values for id.Y
R = calculateR(valX1, valX2, 3, 1); % Calculate R matrix for validation data, for m=1



valY = R * theta;  % validation Yh
MSE_min = sum((Yn_val - valY) .^ 2) / length(Yn_val); % minimum MSE for validation data
MSEs_id(MSE_index) = sum((Yn_id - idY) .^ 2) / length(Yn_id); % MSE for identification data
MSEs_val(MSE_index) = MSE_min;
MSE_index = MSE_index + 1;

%% Iterate through diffeerent values of m (m = 2 to 25)

for m = 2 : 25

    R = zeros(length(X1) * length(X2), size(R, 2) + m + 1); % Update th R matrix with additional columns based on the new m
    R = calculateR(X1, X2, size(R, 2), m);
    theta = R \ Yn_id;  % Estimate the new theta
    idY = R * theta;  % Calculate the fitted values for identification data (id.Y)

    % Calculate the validation predictions and MSE
    valR = calculateR(valX1, valX2, size(R, 2), m);
    valY = valR * theta;

    % Update the minimum MSE
    if (sum((Yn_val - valY) .^ 2) / length(Yn_val)) < MSE_min
        MSE_min = sum((Yn_val - valY) .^ 2) / length(Yn_val);
        m_min = m;
        l_min = size(R, 2);
    end

    % Record the MSEs for both identification and validation datasets
    MSEs_id(MSE_index) = sum((Yn_id - idY) .^ 2) / length(Yn_id);
    MSEs_val(MSE_index) = sum((Yn_val - valY) .^ 2) / length(Yn_val);
    MSE_index = MSE_index + 1;

end

%% Plot the Yh approximation to identification data using a mesh plot
clc
figure
subplot(121);
mesh(id.X{1}, id.X{2}, id.Y); % Original identification data mesh plot
grid on
hold
subtitle('Identification data');

% Calculate Yh using the optimal m and l values found
R = calculateR(X1, X2, l_min, m_min);
theta = R \ Yn_id;
Yh = R * theta;

% Convert Yh vector into a matrix for plotting 
Yh_matrix = zeros(41, 41);
k = 1;
for i = 1 : 41
    for j = 1 : 41
        Yh_matrix(i, j) = Yh(k);
        k = k + 1;
    end
end

% Plot the approximation of the identification data
subplot(122);
mesh(id.X{1}, id.X{2}, Yh_matrix);
grid on;
subtitle('Approximation of Identification data');

% Plot the MSEs for both validation and identification data
figure
plot(MSEs_val(4:20), 'r-*');hold all
plot(MSEs_id(4:20), 'g-*');
plot(1, MSEs_id(4), 'b*');
plot(1, MSE_min, 'b*');
grid on;
title("Visualization of the MSEs for identification and validation");
legend("MSEs_{val}","MSEs_{id}","MSE_{min}");


% Initial plot of the MSEs
figure
plot(MSEs_val, 'r-*');hold all
plot(MSEs_id, 'g-*');
title("Initial plot of the MSEs");
legend("MSEs_{val}", "MSEs_{id}");
grid on;
  
%% Plot the mesh for validation data and approximation
clc
figure
subplot(121);
mesh(valX1, valX2, val.Y); % Original validation data meesh plot
grid on
hold
subtitle('Validation data');

% Calculate Yh using optimal m and l values for validation data
R = calculateR(X1, X2, l_min, m_min);
theta = R \ Yn_id;
Rval = calculateR(valX1, valX2, l_min, m_min);
Yh = Rval * theta;

% Convert the Yh vector to a matrix for plotting
Yh_matrix = zeros(71, 71);
k = 1;
for i = 1 : 71
    for j = 1 : 71
        Yh_matrix(i, j) = Yh(k);
        k = k + 1;
    end
end

% Plot th approximation of validation data
subplot(122);
mesh(valX1, valX2, Yh_matrix);
grid on
subtitle('Approximation of Validation data');

%% Function used to calculate the regressors
function R = calculateR(X1, X2, l,  m)
    rows = 1;
    R = zeros(length(X1) * length(X1), l);


    for i = 1 : length(X1)
        for j = 1 : length(X2)
            R(rows, 1) = 1;
            % disp(['Term: 1']);
            k = 2;   

            % Compute the powers of X1 and X2 up to degree m
            for index = 1 : m
                R(rows, k) = X1(i) ^ (index);
                % disp(['Term: X1^', num2str(k)]);
                k = k + 1;
                R(rows, k) = X2(j) ^ (index);
                % disp(['Term: X2^', num2str(k)]);
                k = k + 1;
            end

            % Compute the cross terms for polynomial fitting
            for a = m - 1 : -1 : 1
                for b = 1 : m - 1
                    if a + b <= m
                        R(rows,k) = (X1(i) ^ a) * (X2(j) ^ b);
                        % disp(['Term: (X1^', num2str(a), ') * (X2^', num2str(b),')']);
                        k = k + 1;
                    end 
                end 
            end

            rows = rows + 1;
        end
    end
end

