% Global variables
d = 0.025;
R = 1;
V = 1;

% Section a
fprintf("\n*** Section a ***\n");
section_a_b(V/2, -V/2, d, R/2, R);

% Section b
fprintf("\n*** Section b ***\n");
section_a_b(V/2, -V/2, d, R/5, R);

% Section c
fprintf("\n*** Section c ***\n");
section_c(V/2, -V/2, d, R);

% Section d
fprintf("\n*** Section d ***\n");
section_d(V, d, R/2, R);


function C_theo = C_theo(R, D)
    % Calculate theoretical capacity for Parallel-plate capacitor
    % assuming D<<R
    C_theo = (epsilon_0 * pi * R^2) / (D);
end

function section_c = section_c(V_top, V_bottom, d, R)
    % Initialize list containing different D values that would be used
     D_list = linspace(d/3, 1, 20);
    % D_list = linspace(0.02, 0.06, 20); Values focusing on small D
    % D_list = linspace(0.1, 1, 10); Values focusing big D
    C_list = zeros(1, length(D_list));
    C_theo_list = zeros(1, length(D_list));
    Diff_list = zeros(1, length(D_list));
    
    % Iterate over all D values at the list, and calculatecapacity and
    % theoretical capacity
    for i=1:length(D_list)
        D = D_list(i);
        C_list(i) = calculate_capacity(V_top, V_bottom, d, D, R);
        C_theo_list(i) = C_theo(R,D); 
        Diff_list(i) = RelativeError(C_list(i), C_theo_list(i));
    end

    % Plot graphs needed for this section
    figure;
    subplot(1,2,1)
    hold on;
    xlabel("D (Distance between plates) [m]");
    ylabel("Capacity [F]");
    plot(D_list, C_list);
    plot(D_list, C_theo_list);
    legend("Calculated Capacity", "Theoretical Capacity");
    hold off;
    
    subplot(1,2,2)
    hold on;
    xlabel("D (Distance between plates) [m]")
    ylabel("Relative Error [%]")
    plot(D_list, 100*Diff_list);
    hold off;

    fprintf('   %s     |  %s  |  %s  \n', ["D","Calculated Capacity", "Theoretical Capacity"]');
    fprintf('%f |    %e       |    %e\n', [transpose(D_list),transpose(C_list),transpose(C_theo_list)]');

    fprintf('\n\n     %s   |  %s\n', ["D","Relative Error"]');
    fprintf('%f |   %f%%\n', [transpose(D_list),transpose(100*Diff_list)]');
end

function section_a_b = section_a_b(V_top, V_bottom, d, D, R)
    % Calculate charge on each palte
    [Q_top, Q_bottom] = CalculateCharge2(d, D, V_top, V_bottom, R);
    % Calculate theoretical Q
    Q_theo = C_theo(R,D) * (V_top - V_bottom);
    % Calculate relative error
    relative_error = RelativeError(Q_top, Q_theo);
    
    % Print Q_top, Q_bottom, relative_error
    fprintf('Q top plate: %e\n', Q_top');
    fprintf('Q bottom plate: %e\n', Q_bottom');
    fprintf('Relative error: %f%%\n', relative_error*100');
end

function section_d = section_d(V, d, D, R)
    % Calculate charge on each palte like section a
    [Q_top, Q_bottom] = CalculateCharge2(d, D, V/2, -V/2, R);
    Q_a = Q_top + Q_bottom;

    % Calculate charge on each palte as requiered at this section
    [Q_top, Q_bottom] = CalculateCharge2(d, D, V, 0, R);
    Q_d = Q_top + Q_bottom;

    % Print Q_a and Q_d
    fprintf('Q top plate: %e\n', Q_top');
    fprintf('Q bottom plate: %e\n', Q_bottom');
    fprintf('Total charge section A: %e\n', Q_a');
    fprintf('Total charge section D: %e\n', Q_d');
end

function C = calculate_capacity(V_top, V_bottom, d, D, R)
    % Calculate charge on each plate
    [Q_top, Q_bottom] = CalculateCharge2(d, D, V_top, V_bottom, R);
    
    % Calculate capacity by fromula: C=Q/V
    C = Q_top / (V_top - V_bottom);
end

function relative_error = RelativeError(actual, theo)
    relative_error = (abs(actual - theo)) / (theo);
end

function [Q_top, Q_bottom] = CalculateCharge2(d, D, V_top, V_bottom, R)
    % Create L describing relation between elements on the top plate
    [x_top,y_top] = get_x_y_range(R,d);
    x_y_combination_top = create_x_y_matrix(x_top,y_top);
    x_y_without_zeros_top = remove_items_out_of_circle(x_y_combination_top, R);
    L_AA = create_l_matrix(x_y_without_zeros_top, V_top, d);

    % Create L describing relation between elements on the bottom plate
    [x_bottom,y_bottom] = get_x_y_range(R,d);
    x_y_combination_bottom = create_x_y_matrix(x_bottom,y_bottom);
    x_y_without_zeros_bottom = remove_items_out_of_circle(x_y_combination_bottom, R);
    L_BB = create_l_matrix(x_y_without_zeros_bottom, V_bottom, d);
    
    % Create L describing relation between elements on the top plate
    % and elements on the bottom plate
    L_AB = create_cross_l_matrix(x_y_without_zeros_top, x_y_without_zeros_bottom, D);
    
    % Create L describing relation between elements on the bottom plate
    % and elements on the top plate
    % Note: Because of the symmetry at the problem there is no need to
    % calculate again
    L_BA = transpose(L_AB);

    % Combine all sub-matrices to a single matrix
    L = [L_AA L_AB; L_BA L_BB];
    
    % Calculate charge on each plate 
    [Q_top, Q_bottom] = calculate_q_top_and_bottom(L, x_y_without_zeros_top, x_y_without_zeros_bottom, V_top, V_bottom);
    
end

function [x, y] = get_x_y_range(R, d)
    % Get combination of all needed x,y values inside the matrix
    % Note: Calculation is a bit different between odd and even
    % number of squares
    num_of_squares = (2*R/d)^2;
    limit = floor(R/d);
    x = -limit : +limit;
    x = x * d;
    
    y = -limit : +limit;
    y = y * d;
    if mod(num_of_squares,2) == 0
        % Fix values in case number of squers is even
        x = x + d/2; 
        x(end) = [];
        y = y + d/2;
        y(end) = []; 
    end
end
function x_y_combination = create_x_y_matrix(x,y)
    % Create a list of the centers of all squares (x,y) at the matrix
    x_y_combination = {y,x};
    [x_y_combination{:}]=ndgrid(x_y_combination{:});
    n=length(x_y_combination);
    x_y_combination = reshape(cat(n+1,x_y_combination{:}),[],n);
end

function x_y_without_zeroes = remove_items_out_of_circle(x_y_combination, R)
    % Remove all items which thier center is outside of the circle
    x_y_without_zeroes = [];
    for i = 1:length(x_y_combination)
        distance = (x_y_combination(i,1)^2 + x_y_combination(i,2)^2);
        if distance < R^2
            x_y_without_zeroes(end+1, 1) = x_y_combination(i,1);
            x_y_without_zeroes(end, 2) = x_y_combination(i,2);
        end
    end
end

function l = create_l_matrix(x_y_without_zeros, V, d)
    % Calculate l matrix as specified at the moments tutorial
    l = zeros(length(x_y_without_zeros),length(x_y_without_zeros));
    % Iterate over all elements twice to find to impact of each element
    % over each element
    for n = 1:length(x_y_without_zeros)
        for m = 1:length(x_y_without_zeros)
            if m == n
                % Impact of the element on it self
                l(m,n) = 0.8814 * 4 /d;
            else
                % Impact of 2 different elements on each other
                d_x = (x_y_without_zeros(n,1) - x_y_without_zeros(m,1))^2;
                d_y = (x_y_without_zeros(n,2) - x_y_without_zeros(m,2))^2;
                l(m,n) = 1/(sqrt(d_x + d_y));
            end
        end
    end
end

function L_AB = create_cross_l_matrix(x_y_without_zeros_top, x_y_without_zeros_bottom, D)
    % Calculate L Matrix between 2 different plates
    L_AB = zeros(length(x_y_without_zeros_bottom),length(x_y_without_zeros_bottom));
    
    % Iterate over all elements twice to find to impact of each element
    % over each element
    for n = 1:length(x_y_without_zeros_top)
        for m = 1:length(x_y_without_zeros_bottom)
            d_x = (x_y_without_zeros_top(n,1) - x_y_without_zeros_bottom(m,1))^2;
            d_y = (x_y_without_zeros_top(n,2) - x_y_without_zeros_bottom(m,2))^2;
            d_z = (D)^2;
            L_AB(n,m) = 1/(sqrt(d_x + d_y + d_z));
        end
    end
end

function [Q_top, Q_bottom] = calculate_q_top_and_bottom(L, x_y_without_zeros_top, x_y_without_zeros_bottom, V_top, V_bottom)
    % Create vector of potentials for each plate and combine them
    V_vec_top = transpose(ones(1, length(x_y_without_zeros_top)) * V_top);
    V_vec_bottom = transpose(ones(1, length(x_y_without_zeros_bottom)) * V_bottom);
    V_vec = [V_vec_top;V_vec_bottom];

    % Solve matrix and calculate charge for each plate
    x_j = linsolve(L, V_vec);
    Q_j = 4 * pi * epsilon_0 * x_j;
    Q_top = sum(Q_j(1:length(Q_j)/2));
    Q_bottom = sum(Q_j(((length(Q_j)/2) +1):length(Q_j)));
end

function epsilon_0 = epsilon_0
    % Note: Unites are [m]
    epsilon_0 = 8.85418*10^-12;
end