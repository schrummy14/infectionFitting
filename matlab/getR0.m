function getR0(data)
    
    R0 = zeros(length(data(:,2))-2,2);
    for k = 2:3
        aa = data(:,k);
        aa = aa(2:end)./aa(1:end-1);
        aa = aa(2:end)./aa(1:end-1);
        R0(:,k-1) = aa;
    end
    
    r0 = R0(:,1);
    d0 = R0(:,2);
    
    plot(data(3:end,1),[r0,d0])
    
    legend('infected R0', 'death R0')
    
end