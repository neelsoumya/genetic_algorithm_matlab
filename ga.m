function [P,best] = ga(pop_size,chrom_len,pm,pc,max_gen)
% /---------------------------------------------------------------------\
% | Copyright (C) 2009  George Bezerra                                   |
% |                                                                      |
% | This program is free software: you can redistribute it and/or modify |
% | it under the terms of the GNU General Public License as published by |
% | the Free Software Foundation, either version 3 of the License, or    |
% | (at your option) any later version.                                  |
% |                                                                      |
% | This program is distributed in the hope that it will be useful,      |
% | but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
% | GNU General Public License for more details.                         |
% |                                                                      |
% | You should have received a copy of the GNU General Public License    |
% | along with this program.  If not, see <http://www.gnu.org/licenses/>.|
% \ ---------------------------------------------------------------------/ 
%
% Inputs:
% pop_size  => population size
% chrom_len => chromosome length
% pm        => probability of mutation
% pc        => probability of crossover
% max_gen   => maximum number of generations
%
% Outputs:
% P    => population
% best => best individual of the population
%
% suggested run: [P,best] = ga(100,100,0.01,0.5,200);


% INITIALIZE POPULATION
P = initialize(pop_size,chrom_len);

% EVALUATION
fit = maxones_fitness(P);

gen = 1;
while gen<=max_gen & max(fit)<chrom_len
    
    % SELECTION
    P = tournament_selection(P,fit,2);
    
    % CROSSOVER
    P = two_point_crossover(P,pc);
    
    % MUTATION
    P = point_mutation(P,pm);
    
    % EVALUATION
    fit = maxones_fitness(P);
   
    % record data
    max_fit(gen) = max(fit);
    mean_fit(gen) = mean(fit);
    
    % display information
    disp_info(gen,max_fit(gen));
    
    gen = gen+1;
end

disp(sprintf('Generation: %d',gen));
disp(sprintf('Best fitness: %d\n',max_fit(end)));

% plot evolution curve
plot(1:length(max_fit), max_fit,'b');
hold on;
plot(1:length(mean_fit), mean_fit,'g');
hold off;
xlabel('Generations');
ylabel('Fitness');
legend('Best fitness','Average fitness','Location','SouthEast');

% output best individual
[m,ind] = max(fit);
best = P(ind,:);


function [P] = initialize(pop_size,chrom_length)
P = round(rand(pop_size,chrom_length));

function [fit] = maxones_fitness(P)
for i=1:size(P,1)
    fit(i) = length(find(P(i,:)));
end

function [P_new] = tournament_selection(P,fit,tourn_size)
for i=1:size(P,1)
    t = ceil(rand(1,tourn_size)*size(P,1));
    [max_fit,winner] = max(fit(t));
    P_new(i,:) = P(t(winner),:);
end

function [P_new] = two_point_crossover(P,pc)
mating_list = randperm(size(P,1));
P_new = [];
while ~isempty(mating_list)
    pair = mating_list(1:2);
    mating_list(1:2) = [];
    if rand<pc
        crossover_points = ceil(rand(1,2)*(size(P,2)));
        point1 = min(crossover_points);
        point2 = max(crossover_points);
        individual1 = P(pair(1),:);
        individual2 = P(pair(2),:);
        individual1(point1:point2) = P(pair(2),point1:point2);
        individual2(point1:point2) = P(pair(1),point1:point2);
        P_new = [P_new;individual1;individual2];
    else
        P_new = [P_new;P(pair,:)];
    end
end

function [P_new] = point_mutation(P,pm)
r = rand(size(P));
mutation_list = find(r<pm);
P_new = P;
P_new(mutation_list(find(P(mutation_list)==1))) = 0;
P_new(mutation_list(find(P(mutation_list)==0))) = 1;

function [] = disp_info(gen,fit)
if mod(gen,10)==0 
    disp(sprintf('Generation: %d',gen));
    disp(sprintf('Best fitness: %d\n',fit));
end

function [P_new] = roulette_selection(P,fit)
fit = (fit - min(fit)).^2;
fit = cumsum(fit);
fit = fit/max(fit);
P_new = [];
for i=1:size(P,1)
    f = find(fit>rand);
    P_new = [P_new;P(f(1),:)];
end
