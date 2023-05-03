%% Laboration 2 Sf2863

%Given constants from assignment
T = [6 8 14 25 12 18 33 8 12];
Lambda = [50 40 45 50 25 48 60 35 15]*1/1000;
c = [12 14 21 20 11 45 75 30 22];
S = [0 0 0 0 0 0 0 0 0];
K=0;

%Start EBO = lambda*T' = 6,1
EBO = [Lambda*T'];

%the cost is equal to the sparepart * cost of the spareparts
Kostnad = [c*S']; 

%runs a while loop until cost exceeds the budget
R_values = linspace(1,9,9);
while Kostnad < 500
    R = [];
    for i = R_values
        R = [R 1 - poisscdf(S(R_values(i)), T(R_values(i))*Lambda(R_values(i)))];
        size(R);
    end
    Comparison = R./c;
    Largest_value_in_column_index = find(Comparison == max(Comparison));
    S(Largest_value_in_column_index) = S(Largest_value_in_column_index) +1;
    S_value = S(Largest_value_in_column_index);
    
    %to obtain a better picture we added the last effecient point after 500
        if c*S' > 506
            break
        else
             Kostnad(end+1) = c*S'; 
             EBO(end+1) = EBO(end) - (1 - poisscdf(S_value-1, ...
                 T(Largest_value_in_column_index)*Lambda(Largest_value_in_column_index)));
        end

 
end

hold on
 plot(Kostnad,EBO,'b-o')
 xlabel('COST')
 ylabel('EBO')
 title('Comparing Marginal allocation with Dynamic programing')

%Here starts the code on dynamic programing

%Given parameters
Lambda = [50 40 45 50 25 48 60 35 15]/1000; %Arrival rates of malfunctioning LRU's
T = [6 8 14 25 12 18 33 8 12];              %Repair time for malfunctioning LRU
C = [12 14 21 20 11 45 75 30 22];           %Purchase cost per LRU

%Creating vectors that will be used
new_values = [];     %Vector that will contain the optimal values
X_n = [];   %Vector that will contain when all decisions are made

current_LRU = 9; %the ninth stage
number_of_spareparts = 0; %amount of spareparts of type n

%Creating the values for the last LRU', when n = 9. 
max_budget = 500;
b = 0;

%the trivial case, running a while loop from 0 to 500. 
while b <= max_budget 
     budget = b; 
     if budget/C(current_LRU) == 0
         number_of_spareparts = 0;
         new_values = [new_values; EBO_function(number_of_spareparts, current_LRU)];
         X_n = [X_n; number_of_spareparts];
     else
         number_of_spareparts = fix(budget/C(current_LRU));
         new_values = [new_values; EBO_function(number_of_spareparts, current_LRU)];
         X_n = [X_n; number_of_spareparts]; %The decision vector
     end
     b = b + 1;
     
end
size(b);

optimal_value_vector = [new_values]; %The optimal value vector will be a 501x9 matrix containing all the information

%For when n = 8, 7, 6, 5, 4, 3, 2, 1, working recursively 
current_LRU = 8;
budget = 0;
number_of_spareparts = 0;
x = []; %The decisions
new_values = [];
generated_values = [];

%to run the different stages(LRU's) recursively we created a while loop
budget_values = linspace(0, max_budget, 501);
%start to check if we runned through all of the LRU's 
while current_LRU > 0
    for b = budget_values
        budget = budget_values(b+1);
        %stop condition if we exceeded budget
        while budget >= number_of_spareparts*C(current_LRU)
            generated_values = [generated_values EBO_function(number_of_spareparts, current_LRU) + optimal_value_vector(budget-number_of_spareparts*C(current_LRU)+1, 1)];
            number_of_spareparts = number_of_spareparts + 1;
        end
        to_insert = min(generated_values);
        new_values = [new_values; to_insert];
        new_decision = find(generated_values == min(generated_values)) - 1;
        if length(new_decision) > 1
            new_decision = new_decision(1);
        end
        x = [x; new_decision];
        number_of_spareparts = 0;
        generated_values = [];
    end
    optimal_value_vector = [new_values optimal_value_vector];
    X_n = [x X_n];
    new_values = [];
    x = [];
    X_n;
    current_LRU = current_LRU - 1;
    
end


%the budget below was given by assignment
assignment_budget = [0 100 150 350 500]';
budget = 500;
cost = 0;
total_costs = [];
strategies = [];
optimal_path = [];

EBO = [];
current_index = 1;
%need to run this 5 times for the values of 0 100 150 350 500
while current_index <= length(assignment_budget)
    %setting budget to the budget we got from the question
    budget = assignment_budget(current_index);
    EBO = [EBO; optimal_value_vector(budget+1,1)];
    for b = 1:9
        current_LRU = b;
        allocation = X_n(budget+1,current_LRU);
        cost = cost + allocation*C(current_LRU);
        strategies = [strategies allocation];
        budget = budget-allocation*C(current_LRU);
    end
    current_index = current_index + 1;
optimal_path = [optimal_path; strategies];
total_costs = [total_costs; cost];
strategies = [];
cost = 0;

end

%to be able to se the table, take away the percentage sign "%"
table(optimal_path(:,1), optimal_path(:,2), optimal_path(:,3), optimal_path(:,4), optimal_path(:,5), optimal_path(:,6), optimal_path(:,7), optimal_path(:,8), optimal_path(:,9), total_costs, EBO, assignment_budget)
hold on
plot(total_costs,EBO, 'Color', 'red')
plot(total_costs,EBO,'*', 'LineWidth', 4,  'Color', 'Red')
legend({'y = Marginal Allocation','y = Dynamic Programming'},'Location','southwest')
hold off










