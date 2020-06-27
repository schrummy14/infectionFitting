function compStates(s1, s2, pop1, pop2, stateName1, stateName2)

    if nargin < 6
        stateName1 = 'State 1';
        stateName2 = 'State 2';
        if nargin < 4
            pop1 = 1;
            pop2 = 1;
        else
            pop1 = pop1/100;
            pop2 = pop2/100;
        end
    else
        pop1 = pop1/100;
        pop2 = pop2/100;
    end
    x1 = s1(:,1) - s1(1,1);
    i1 = s1(:,2)/pop1;
    d1 = s1(:,3)/pop1;
    
    x2 = s2(:,1) - s2(1,1);
    i2 = s2(:,2)/pop2;
    d2 = s2(:,3)/pop2;
    
    figure(2)
    subplot(2,1,1)
    plot(x1,i1,'b-',x2,i2,'r--')
%     semilogy(x1,i1,'b-',x2,i2,'r--')
    if nargin < 4
        ylabel('Number of Infections')
    else
        ylabel('Percent Infected')
    end
    legend(stateName1,stateName2,'Location','NW')
    subplot(2,1,2)
    plot(x1,d1,'b-',x2,d2,'r--')
%     semilogy(x1,d1,'b-',x2,d2,'r--')
    xlabel('Days')
    if nargin < 4
        ylabel('Number of Deaths')
    else
        ylabel('Percent Deaths')
    end
    legend(stateName1,stateName2,'Location','NW')
    
end
    
    