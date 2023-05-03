%%       SF2863 Home assignment 1, Alexander Råberg and Jakob Amaya Scott        


%Note : All code for the assignment is in this file. The file has been
%divided into different subections. by running the first section of the
%code, the result(s) from question 1 will be obtained, and by running
%section 2, the result(s) from question 2 will be obtained, etc. 


%%        Question 2: Determine the intensity matrix

format short

%Given parameters
v = 18;
v1 = 15;
v2 = 10;
lambda1 = 14;
lambda2 = 20;
my1 = 10;
my2 = 8;
h = 0.001;
n = 3; %Number of employees

%Strategy 1, if 1 engine is broken, then all 3 workers will be repairing
%the broken machine. If both engines are broken, then all three workers will
%be repairing machine 1

%Creating the intensity matrix
Q = zeros(4,4);
%First row
Q(1,1) = -(lambda1 + lambda2);
Q(1,2) = lambda1;
Q(1,3) = lambda2;
Q(1,4) = 0;
%Second row
Q(2,1) = my1*n;
Q(2,2) = -(my1*n+lambda2);
Q(2,3) = 0; 
Q(2,4) = lambda2;
%Third row
Q(3,1) = n*my2;
Q(3,2) = 0;
Q(3,3) = -(n*my2 + lambda1);
Q(3,4) = lambda1;
%Fourth row. This row in the matrix will change for the different
%strategies. For this strategy (strategy 1), the third element in the
%fourth row will be = 30, which is the intensity if all 3 workers are
%working on machine 1. 
Q(4,1) = 0;
Q(4,2) = 0;
Q(4,3) = n*my1;
Q(4,4) = -n*my1;

Q_1 = Q;
intensity_matrix_strategy_1 = Q 

%Strategy 2, if 1 engine is broken, then all 3 workers will be repairing
%the broken machine. If both engines are broken, 2 workers will
%be repairing machine 1 and 1 worker will be repairing machine 2
Q_2 = Q; 
Q_2(4, 1:end) = [0 8 20 -28];
intensity_matrix_strategy_2 = Q_2

%Strategy 3, if 1 engine is broken, then all 3 workers will be repairing
%the broken machine. If both engines are broken, 1 worker will
%be repairing machine 1 and 2 workers will be repairing machine 2
Q_3 = Q; 
Q_3(4, 1:end) = [0 16 10 -26];
intensity_matrix_strategy_3 = Q_3

%Strategy 4, if 1 engine is broken, then all 3 workers will be repairing
%the broken machine. If both engines are broken, all three workers will be
%repairing machine 2
Q_4 = Q; 
Q_4(4, 1:end) = [0 24 0 -24];
intensity_matrix_strategy_4 = Q_4


%% Question 4 & 5, determine the stationary distributions and their corresponding average speed
format long
%Stationary distribution for strategy 1
Q_1_transpose = Q_1';       %Transposing the matrix so that we can solve for the columns.   
Q_1_transpose(4,1:4)=ones;  %The condition that sum of all pi's need to = 1

to_solve = [0;0;0;1];
stationary_distribution = Q_1_transpose\to_solve; %Our solution, i.e. the stationary distribution
speed = stationary_distribution.*[v, v1,v2,0]';   %Speed for first strategy
average_speed = sum(speed);                       %The average speed

stationary_distribution = stationary_distribution';

stationary_distribution_1 = stationary_distribution
speed_1 = speed;
average_speed_1 = average_speed

%Stationary distribution for strategy 2
Q_2_transpose = Q_2';       %Transposing the matrix so that we can solve for the columns.   
Q_2_transpose(4,1:4)=ones;  %The condition that sum of all pi's need to = 1

to_solve = [0;0;0;1];
stationary_distribution_2 = Q_2_transpose\to_solve; %Our solution, i.e. the stationary distribution
speed_2 = stationary_distribution_2.*[v, v1,v2,0]'; %Speed for second strategy
average_speed_2 = sum(speed_2);                     %The average speed

stationary_distribution_2 = stationary_distribution_2'
speed_2;
average_speed_2

%Stationary distribution for strategy 3
Q_3_transpose = Q_3';       %Transposing the matrix so that we can solve for the columns.   
Q_3_transpose(4,1:4)=ones;  %The condition that sum of all pi's need to = 1

to_solve = [0;0;0;1];
stationary_distribution_3 = Q_3_transpose\to_solve; %Our solution, i.e. the stationary distribution
speed_3 = stationary_distribution_3.*[v, v1,v2,0]'; %Speed for third strategy
average_speed_3 = sum(speed_3);                     %The average speed

stationary_distribution_3 = stationary_distribution_3'
speed_3;
average_speed_3

%Stationary distribution for strategy 4
Q_4_transpose = Q_4';       %Transposing the matrix so that we can solve for the columns.   
Q_4_transpose(4,1:4)=ones;  %The condition that sum of all pi's need to = 1

to_solve = [0;0;0;1];
stationary_distribution_4 = Q_4_transpose\to_solve; %Our solution, i.e. the stationary distribution
speed_4 = stationary_distribution_4.*[v, v1,v2,0]'; %Speed for fourth strategy
average_speed_4 = sum(speed_4);                     %The average speed

stationary_distribution_4 = stationary_distribution_4'
speed_4;
average_speed_4

%% Question 6a) & 7a): How we determine the time to the next jump and where to jump next


%Creating the P-matrix that will be used to determined where to jump next.
%For the strategy 1
probability_transition_matrix = zeros(4, 4);
for i = 1:4    
    current_row = Q(i, 1:4); 
    index = find((current_row) > 0);  %Want to find all elements in the row greater than 0
    probability_transition_matrix(i, index) = current_row(index)./(-Q(i, i)); %Dividing all non-negative elements by Q(i, i) to obtain the transition probability
    positive = find(index > 0);    
end
probability_transition_matrix


%% Question 6a) & 7a): How we determine the time to the next jump and where to jump next

T = 1000;            %The time that we will run the simulation, can the T to what we want
starting_state = 1; %Start the simulation in state 1, meaning that both engines are working
current_state = starting_state;
total_time = 0; %How long time the simulation has run for

time_in_each_state = zeros(1, 4); %Matrix that will contain all the time spent in each state
number_jumps = 1;            %How many

while total_time <= T 
    index = Q(current_state, 1:4);
    positive = find(index > 0);            
     
    %How we determine our next jump, using exprnd. The input parameter will
    %be 1/(-Q(current_state, current_state), meaning the diagonal element
    %in the Q-matrix
    time_interval = exprnd(1/-(Q(current_state, current_state))); 
    
    %Determining the next jump by using the transition probabilities in the
    %transition probability matrix above. We return the index of the next
    %state that we jump to. 
    next_state = randsrc(1,1,[1,2, 3, 4;probability_transition_matrix(current_state, 1:4)]);    
    
    time_in_each_state(number_jumps, current_state) = time_interval; %Updating the vector with values
    total_time = total_time + time_interval;           %Updating the total time
    current_state = next_state;
    number_jumps = number_jumps + 1;
end

%Vector used to store the total time that the simulation has spend in each state
total_time_in_each_state = zeros(1, 4); 
for i = 1:4
    total_time_in_each_state(i) = sum(time_in_each_state(:, i));    %Obtaining the total time in each state by adding the columns in the vector
end


simulated_stationary_continuous = total_time_in_each_state./total_time; %How much time we spend in each state as a perfect, i.e
simulated_stationary_continuous

%% Question 8a) Take average to determine estimates of the expected speed of the ferry
average_speed_continuous = sum(simulated_stationary_continuous*[v, v1,v2,0]');
average_speed_continuous


%% Question 6b): Determine the transition matrix of a discrete time Markov chain that will approximate the continuous time process

discrete_Q = zeros(4, 4);   %Creating the matrix for the discrete time approximation of the continuous Markov process

for i = 1:4   
    discrete_Q(i, i) = 1 + h*Q(i, i);   %Want all diagonal elements of matrix to be = 1 + h*Q(i, i)
end

x = [2 3 4 5 7 8 9 10 12 13 14 15];     %Index of all non-diagonal elements of matrix
length(x);
for i = 1:length(x)        
    discrete_Q(x(i)) = discrete_Q(x(i)) + h*Q(x(i));    %Want all non-diagonal elements to be = h*Q(i, i)
end

discrete_Q

%% Question 7b): Use the discretization to simulate the process 

T = 100;            %The time that we will run the simulation, can the T to what we want
starting_state = 1; %Start the simulation in state 1, meaning that both engines are working
current_state = starting_state;
total_time = 0; %How long time the simulation has run for

number_iterations = T/h;
time_interval = h;  %Can move this outside the loop since the time is the for every iteration is the same

number_of_runs = 4; %How many times we will run the program and take averages
average_stationary_distribution_discrete_vector = zeros(number_of_runs, 4);
average_speed_discrete_vector = zeros(number_of_runs, 1);

for j = 1:number_of_runs
    total_time = 0;
    time_in_each_state = zeros(T/h, 4); %Matrix that will contain all the time spent in each state
    number_jumps = 1;            %How many jumps that have been made so far in the simulation
    total_time_in_each_state = zeros(1, 4);     %The percentage of time spent in each state, i.e. the simulated discrete stationary distribution.

    for i=1:number_iterations
        index = discrete_Q(current_state, 1:4);
        positive = find(index > 0);            

        %Determining the next jump by using the transition probabilities in the
        %transition probability matrix above. We return the index of the next
        %state that we jump to. 
        next_state = randsrc(1,1,[1,2, 3, 4; discrete_Q(current_state, 1:4)]);    

        time_in_each_state(number_jumps, current_state) = time_interval; %Updating the vector with values
        total_time = total_time + time_interval;           %Updating the total time
        current_state = next_state;
        number_jumps = number_jumps + 1;
        if i == number_iterations
            for k = 1:4
                total_time_in_each_state(k) = sum(time_in_each_state(:, k))./total_time;    %Obtaining the total time in each state by adding the columns in the vector
            end
            average_stationary_distribution_discrete_vector(j, 1:end) = total_time_in_each_state;

        end
    end
end

estimated_speed_discrete_vector = zeros(length(average_stationary_distribution_discrete_vector(:, 1)), 1);

for i = 1:length(average_stationary_distribution_discrete_vector(:, 1))
    estimated_speed_discrete_vector(i) = sum(average_stationary_distribution_discrete_vector(i, 1:end).*[v v1 v2 0]);
end

error = std(estimated_speed_discrete_vector)

%%  8b): Take averages to determine estimates of the expected speed of the ferry
formatSpec = ('The average speed for %d trials with a T = %4d was %g');
sprintf(formatSpec, number_of_runs, T, mean(estimated_speed_discrete_vector))
%average_speed_discrete = mean(estimated_speed_discrete_vector)


%%  9b): How the computation times between the three different methods compare
format longG

%For the analytic solution

%Stationary distribution for strategy 1
tic
Q_1_transpose = Q_1';       %Transposing the matrix so that we can solve for the columns.   
Q_1_transpose(4,1:4)=ones;  %The condition that sum of all pi's need to = 1

to_solve = [0;0;0;1];
stationary_distribution = Q_1_transpose\to_solve; %Our solution, i.e. the stationary distribution
speed = stationary_distribution.*[v, v1,v2,0]';   %Speed for first strategy
average_speed = sum(speed);                       %The average speed

stationary_distribution = stationary_distribution';

stationary_distribution_1 = stationary_distribution;
speed_1 = speed;
average_speed_1 = average_speed;

analytic_time = toc;

%For the Continuous time approach
T = 100;            %The time that we will run the simulation, can the T to what we want
starting_state = 1; %Start the simulation in state 1, meaning that both engines are working
current_state = starting_state;
total_time = 0; %How long time the simulation has run for

time_in_each_state = zeros(1, 4); %Matrix that will contain all the time spent in each state
number_jumps = 1;            %How many

tic
while total_time <= T 
    index = Q(current_state, 1:4);
    positive = find(index > 0);            
     
    %How we determine our next jump, using exprnd. The input parameter will
    %be 1/(-Q(current_state, current_state), meaning the diagonal element
    %in the Q-matrix
    time_interval = exprnd(1/-(Q(current_state, current_state))); 
    
    %Determining the next jump by using the transition probabilities in the
    %transition probability matrix above. We return the index of the next
    %state that we jump to. 
    next_state = randsrc(1,1,[1,2, 3, 4;probability_transition_matrix(current_state, 1:4)]);    
    
    time_in_each_state(number_jumps, current_state) = time_interval; %Updating the vector with values
    total_time = total_time + time_interval;           %Updating the total time
    current_state = next_state;
    number_jumps = number_jumps + 1;
end
continuous_time = toc;

%For the discretization approach
current_state = starting_state;
total_time = 0; %How long time the simulation has run for


time_in_each_state = zeros(1, 4); %Matrix that will contain all the time spent in each state
number_jumps = 1;            %How many
number_iterations = T/h;

tic
for i = 1:number_iterations 
    index = Q(current_state, 1:4);
    positive = find(index > 0);            
    time_interval = h; 
    
    %Determining the next jump by using the transition probabilities in the
    %transition probability matrix above. We return the index of the next
    %state that we jump to. 
    next_state = randsrc(1,1,[1,2, 3, 4;probability_transition_matrix(current_state, 1:4)]);    
    
    time_in_each_state(number_jumps, current_state) = time_interval; %Updating the vector with values
    total_time = total_time + time_interval;           %Updating the total time
    current_state = next_state;
    number_jumps = number_jumps + 1;
end
discrete_time = toc;

table_size = [1 4];
varTypes = ["double", "double","double","double"];
varNames = ["T", "Analytic time","Continuous time","Discrete time"];
information = [T analytic_time continuous_time discrete_time];

temps = table('Size',table_size,'VariableTypes',varTypes,'VariableNames',varNames);
temps(1, :) = {information(1), information(2), information(3), information(4)}


%%  Question 10:        Average time to full breakdown, analytical

%Creating the matrix that will be used to solve the equation system for the
%average time until breakdown. This will give us the analytic time. 
t_matrix = [-1 Q(1, 2)/abs(Q(1, 1)) Q(1, 3)/abs(Q(1, 1));
            Q(2, 1)/abs(Q(2, 2)) -1 Q(2, 3)/abs(Q(2, 2));
            Q(3,1)/abs(Q(3, 3)) Q(3, 2)/abs(Q(3, 3)) -1];

to_solve = [-1/abs(Q(1, 1)) -1/abs(Q(2, 2)) -1/abs(Q(3, 3))]';

t_vector = t_matrix\to_solve;  %This will give us [t0, t1, t2] for the analytic solution

analytic_time_until_breakdown = t_vector(1);

format longG
number_of_attempts = [10 100 1000 10000 100000];  %For the different number of iterations that we will compare

%total_simulated_time_vector = zeros(1, number_of_attempts);
%number_of_steps_vector = zeros(1, length(number_of_attempts));

mean_time_until_breakdown_vector = zeros(length(number_of_attempts), 1);



for k = 1:length(number_of_attempts)
    
    
%    number_of_steps_vector = zeros(1, number_of_attempts(k));
    total_simulated_time_vector = zeros(1, number_of_attempts(k));
    length(total_simulated_time_vector);
 
    
    for j = 1:number_of_attempts(k)
        
        current_state = starting_state;
        number_jumps = 1;
        loop = true;
        total_time = 0;
        time_in_each_state = zeros(1, 4);

            while loop == true             

                time_interval = exprnd(1/-(Q(current_state, current_state))); 
                next_state = randsrc(1,1,[1, 2, 3, 4;probability_transition_matrix(current_state, 1:4)]);

                time_in_each_state(number_jumps, current_state) = time_interval; %Updating the vector with values
                total_time = total_time + time_interval;

                if next_state == 4
                   total_simulated_time_vector(j) = total_time;
                   number_of_steps_vector(j) = number_jumps;
                   mean_time_until_breakdown_vector(k) = mean(total_simulated_time_vector);
                   loop = false;
                   

                end
                current_state = next_state;
                number_jumps = number_jumps + 1;
            end
    end
end

mean_time_until_breakdown_vector;


table_size = [length(mean_time_until_breakdown_vector) 3];
varTypes = ["double","double", "double"];
varNames = ["Number of runs", "Average time until first breakdown", "Difference from theoretical value"];

%information = ["10 100"]
rowNames =   ["10", "100", "1 000", "10 000", "100 000"]';

temps = table('Size',table_size,'VariableTypes',varTypes,'VariableNames',varNames, ...
    'RowNames', rowNames);
diff_vector = abs(analytic_time_until_breakdown - mean_time_until_breakdown_vector);

temps(1, :) = {number_of_attempts(1), mean_time_until_breakdown_vector(1), diff_vector(1)};
temps(2, :) = {number_of_attempts(2), mean_time_until_breakdown_vector(2), diff_vector(2)};
temps(3, :) = {number_of_attempts(3), mean_time_until_breakdown_vector(3), diff_vector(3)};
temps(4, :) = {number_of_attempts(4), mean_time_until_breakdown_vector(4), diff_vector(4)};
temps(5, :) = {number_of_attempts(5), mean_time_until_breakdown_vector(5), diff_vector(5)}