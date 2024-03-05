% Global variables
V = 1;
R = 1;
d = [0.25,0.15,0.12,0.1,0.075,0.05,0.025,0.02];
Q = zeros(1, length(d));
Q_theo = 8*epsilon_0 * V * R;

% Iterate over all d values and calculate total charge
for i=1:length(d)
    Q(i) = CalculateCharge1(d(i), V, R);
end

% Plot graph mapping between d and total charge 
figure;
hold on
xlabel("d (length of each square) [m]");
ylabel("Total charge [C]");
ylim([6.8*10^-11 7.2*10^-11]);
yline(Q_theo);
plot(d, Q);
legend("Theoretical charge from section A", "Charge calculated by Moments method");
hold off

fprintf('Theoretical charge from section A: %e\n', Q_theo');
fprintf('     %s   |  %s\n', ["d","Q"]');
fprintf('%f|%e\n', [transpose(d),transpose(Q)]');



function Q = CalculateCharge1(d, V, R)
    % Create list of all squares centers, and remove items outside of
    % circle
    [x,y] = get_x_y_range(R,d);
    x_y_combination = create_x_y_matrix(x,y);
    x_y_without_zeros = remove_items_out_of_circle(x_y_combination, R);
    
    % Create l matrix as specified at the moments tutorial
    l = create_l_matrix(x_y_without_zeros, V, d);
    
    % Calculate total charge
    Q = calculate_q(l, x_y_without_zeros, V);
end

function [x, y] = get_x_y_range(R, d)
    % Get combination of all needed x,y values inside the matrix
    % Note: Calculation is a bit different between odd and even
    % number of squares
    num_of_squers = (2*R/d)^2;
    limit = floor(R/d);
    x = -limit : +limit;
    x = x * d;
    
    y = -limit : +limit;
    y = y * d;
    if mod(num_of_squers,2) == 0
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

function Q = calculate_q(l, x_y_without_zeros, V)
    % Create vector full of V values
    V_vec = transpose(ones(1, length(x_y_without_zeros)) * V);     
    % Solve matrix to find x_j values
    x_j = linsolve(l, V_vec);
    % Calculate charge of each element
    Q_j = 4 * pi * epsilon_0 * x_j;
    % Sum all elements to receive total charge
    Q = sum(Q_j);
end

function epsilon_0 = epsilon_0
    % Note: Unites are [m]
    epsilon_0 = 8.85418*10^-12;
end