% Solution of fractional-order delay logistic equation based Matlab software

% References: please see the references[5,7-9,10] in the  paper "Chaos
% detection and parameter identification in fractional-order chaotic systems with delay"

%Email: liguoychina@gmail.com



% clear all

clc

close all

h=0.01;

N=12000;

q1=0.9; %the order of modified logistic model

a=26;

r=-53;%parameter coffecients of the equation

Ndelay=50;

x1 = zeros(1,Ndelay+N+1);

x = zeros(1,Ndelay+N+1);

x(1:Ndelay) = 0.5;

x0=x(Ndelay);

x1(Ndelay+1)=x0+h^q1*(-a*x0+r*x(1)-r*x(1)*x(1))/(gamma(q1)*q1);

x(Ndelay+1)=x0+h^q1*(-a*x(Ndelay+1)+r*x(2)-r*x(2)*x(2)+q1*(-a*x0+r*x(1)-r*x(1)*x(1)))/gamma(q1+2);  

wb = waitbar(0.0,'Solving - 0%');
for n=1:N

    M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(-a*x0+r*x(1)-r*x(1)*x(1));

    N1=((n+1)^q1-n^q1)*(-a*x0+r*x(1)-r*x(1)*x(1));

    for j=1:n

        M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(-a*x(Ndelay+j)+r*x(j)-r*x(j)*x(j));

        N1=N1+((n-j+1)^q1-(n-j)^q1)*(-a*x(Ndelay+j)+r*x(j)-r*x(j)*x(j));  

    end

    x1(Ndelay+n+1)=x0+h^q1*N1/(gamma(q1)*q1);

    x(Ndelay+n+1)=x0+h^q1*(-a*x1(Ndelay+n+1)+r*x(j+1)-r*x(j+1)*x(j+1)+M1)/gamma(q1+2);
    
    if mod(n,100)==1
        nOnN = n/N;
        waitbar(nOnN,wb,['Solving - ', num2str(round(100*nOnN)), '%']);
    end

end
delete(wb);

xresult = zeros(1,N-Ndelay+1);

yresult = zeros(1,N-Ndelay+1);

for n=2*Ndelay+1:N+Ndelay+1

   xresult(n-2*Ndelay)=x(n);

   yresult(n-2*Ndelay)=x(n-Ndelay);

end

plot(xresult(500:end),yresult(500:end),'--')

axis square

xlabel('x(t)')

ylabel('x(t-0.5)')