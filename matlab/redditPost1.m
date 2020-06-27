clc
clearvars -except dataIA

% Days from May 1st
may1 = datetime(2020,5,1);
numDays = floor(daysact(may1,datetime()));

% n_day_projection(dataIA(end-10-numDays:end-numDays,:),5,dataIA(end-numDays+1:end,:));

% pause

n_day_projection(dataIA(end-10:end,:),5)