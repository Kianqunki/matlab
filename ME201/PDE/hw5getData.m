file = 'Prairie Grass.xls';

% read table 2 from the excel file and save it
table2 = xlsread(file);
save data table2;