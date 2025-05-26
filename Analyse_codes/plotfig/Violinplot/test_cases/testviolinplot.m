
% TEST CASE 1
disp('Test 1: Violin plot default options');
load carbig MPG Origin
Origin = cellstr(Origin);
figure
vs = violinplot(MPG, Origin);
ylabel('Fuel Economy in MPG');
xlim([0.5, 7.5]);

disp('Test 1 passed ok');

% TEST CASE 2
disp('Test 2: Test the plot ordering option');
grouporder={'USA','Sweden','Japan','Italy','Germany','France','England'};
    
figure;
vs2 = violinplot(MPG,Origin,'GroupOrder',grouporder);
disp('Test 2 passed ok');

%other test cases could be added here
load carbig MPG Origin % load some example data
Origin = string(Origin);
figure
vp = violinplot(MPG, Origin);
vp(2).ViolinColor = [1 0 0]; % make the second violin red
