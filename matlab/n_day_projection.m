function [yyi5, CIi5, yyd5, CId5] = n_day_projection(data,dayProj,data2)
    
    if nargin < 2
        dayProj = 5;
    end
    x = data(:,1);
    dataInfected = data(:,2);
    dataDeaths = data(:,3);
    x = x - x(1) + 1;
    AlphaVal = 0.01;
    Prediction = 'observation';
    Simultaneous = 'on';
    
    if nargin == 3
        dataInfected2 = data2(:,2);
        dataDeaths2 = data2(:,3);
        x2 = data2(:,1)-data(end,1);
    end
    
%     Prediction = 'curve';
%     Simultaneous = 'off';
    
    n = x; % S(1:length(data))';
    LMi = stepwiselm(x,log(dataInfected),'poly2');
    LMd = stepwiselm(x,log(dataDeaths),'poly2');
    display(LMi);
    display(LMd);
    figure(1)
    clf
    subplot(2,2,1)
    plot(LMi)
    xlabel('Days')
    ylabel('Log(Number Infected)')
    title('Raw Data Infected')
%     figure(2)
    subplot(2,2,2)
    plot(LMd)
    xlabel('Days')
    ylabel('Log(Number of Deaths)')
    title('Raw Data Deaths')
%     figure(3)
    subplot(2,2,3)

    nn = round((1:0.1:n(end)+dayProj)',1);
    [yy, CI] = LMi.predict(nn,...
        'Alpha',AlphaVal, ...
        'Prediction',Prediction, ...
        'Simultaneous', Simultaneous);
    n5 = (n(end)+1:n(end)+dayProj)';
    [yy5,CI5] = (LMi.predict(n5,...
        'Alpha',AlphaVal, ...
        'Prediction',Prediction, ...
        'Simultaneous', Simultaneous));
    
    ids = zeros(dayProj,1);
    ids(1) = find(nn==n5(1));
    for k = 2:dayProj
        ids(k) = ids(k-1) + 10;
    end
    
    disp('CI Low, Prediction, CI High')
    disp(floor(exp([CI5(:,1), yy5, CI5(:,2)])))

    if nargin == 3
        h = plot( ...
            n-length(dataDeaths),dataInfected/1000,'k*', ...
            nn-length(dataDeaths),exp(yy)/1000,'b-o', ...
            x2,dataInfected2/1000,'ks', ...
            nn-length(dataDeaths),exp(CI(:,1))/1000,'r--', ...
            nn-length(dataDeaths),exp(CI(:,2))/1000,'r--');
    else
        h = plot( ...
            n-length(dataDeaths),dataInfected/1000,'k*', ...
            nn-length(dataDeaths),exp(yy)/1000,'b-o',...
            nn-length(dataDeaths),exp(CI(:,1))/1000,'r--', ...
            nn-length(dataDeaths),exp(CI(:,2))/1000,'r--');
    end
    
    for k = 1:length(h)
        h(k).LineWidth = 2;
        h(k).MarkerSize = 8;
    end
    
    h(2).MarkerIndices = ids;
    
    xlabel 'Days from Today'
    ylabel 'Number of Infected (Per Thousand)'
    grid on
    if nargin == 3
        lg = legend( ...
            'Provided Data', ...
            [num2str(dayProj),'-Day-Forecast'], ...
            'Reported Values', ...
            [num2str(100*(1-AlphaVal)),'% CI'], ...
            'Location','NW');
    else
        lg = legend( ...
            'Provided Data', ...
            [num2str(dayProj),'-Day-Forecast'], ...
            [num2str(100*(1-AlphaVal)),'% CI'], ...
            'Location','NW');
    end
    
    yyi5 = floor(exp(yy5));
    CIi5 = floor(exp(CI5));
    
    
    
    subplot(2,2,4)
    [yy, CI] = LMd.predict(nn,...
        'Alpha',AlphaVal, ...
        'Prediction',Prediction, ...
        'Simultaneous', Simultaneous);
    n5 = (n(end)+1:n(end)+dayProj)';
    [yy5,CI5] = (LMd.predict(n5,...
        'Alpha',AlphaVal, ...
        'Prediction',Prediction, ...
        'Simultaneous', Simultaneous));
    
    [maxYY,maxYYid] = max(yy);
    [~,maxYY5id] = max(yy5);
    [maxCI1,maxCI1id] = max(CI(:,1));
    [maxCI2,maxCI2id] = max(CI(:,2));
    
    [~,maxCI51id] = max(CI5(:,1));
    [~,maxCI52id] = max(CI5(:,2));
    
    ids = zeros(dayProj,1);
    ids(1) = find(nn==n5(1));
    for k = 2:dayProj
        ids(k) = ids(k-1) + 10;
    end
    
    yy(maxYYid:end) = maxYY;
    yy5(maxYY5id:end) = maxYY;
    CI(maxCI1id:end,1) = maxCI1;
    CI(maxCI2id:end,2) = maxCI2;
    
    CI5(maxCI51id:end,1) = maxCI1;
    CI5(maxCI52id:end,2) = maxCI2;
    
    disp('CI Low, Prediction, CI High')
    disp(floor(exp([CI5(:,1), yy5, CI5(:,2)])))

    if nargin == 3
        h = plot( ...
            n-length(dataDeaths),dataDeaths/1000,'k*', ...
            nn-length(dataDeaths),exp(yy)/1000,'b-o',...
            x2,dataDeaths2/1000,'ks', ...
            nn-length(dataDeaths),exp(CI(:,1))/1000,'r--', ...
            nn-length(dataDeaths),exp(CI(:,2))/1000,'r--');
    else
        h = plot( ...
            n-length(dataDeaths),dataDeaths/1000,'k*', ...
            nn-length(dataDeaths),exp(yy)/1000,'b-o',...
            nn-length(dataDeaths),exp(CI(:,1))/1000,'r--', ...
            nn-length(dataDeaths),exp(CI(:,2))/1000,'r--');
    end
    
    for k = 1:length(h)
        h(k).LineWidth = 2;
        h(k).MarkerSize = 8;
    end
    
    h(2).MarkerIndices = ids;
    
    xlabel 'Days from Today'
    ylabel 'Number of Deaths (Per Thousand)'
%     title 'Deaths'
    grid on
    
    if nargin == 3
        lg = legend( ...
            'Provided Data', ...
            [num2str(dayProj),'-Day-Forecast'], ...
            'Reported Values', ...
            [num2str(100*(1-AlphaVal)),'% CI'], ...
            'Location','NW');
    else
        lg = legend( ...
            'Provided Data', ...
            [num2str(dayProj),'-Day-Forecast'], ...
            [num2str(100*(1-AlphaVal)),'% CI'], ...
            'Location','NW');
    end
    
    
    yyd5 = floor(exp(yy5));
    CId5 = floor(exp(CI5));
end
    