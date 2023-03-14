clc;
clear all;

%% Unscented Kalman Filter and smoother for the Heston Model
% inputs:
% parameters = [kappa, theta, eta/sigma, rho];
parameters = [3, 0.044664571, 0.019238716, -0.306167243];
% and the stock drift
mu = 0.053;
% plus a vector of prices, Stock
data = readtable("S&P500.csv");
Stock = data.StockPrice;
dt = 1/252;
% mu = stock drift
% kappa = reversion speed
% theta = long run variance
% sigma = volatility of sqrt(var) in varianace equation
% rho = correlation

kappa = parameters(1);
theta = parameters(2);
sigma = parameters(3);
rho = parameters(4);

n = length(Stock);   %Length of stock-price vector

L = 3;              %Dimension

%Default values
alpha = 0.001;      %Default
k = 3-L;            %Deafult
beta = 2;           %Optimal for Gaussian dist
eps = 0.00001;
lambda = alpha^2*(L+k)-L;

%Initialising
X = ones(1,7);
Xa = ones(3,7);
Y = ones(1,7);
estimates = ones(1,n);
u(1) = 0;
v(1) = 1;
vol = 0;

estimates(1) = Stock(1) + eps;
estimates(2) = Stock(1) + eps;

%%Initial parameter guesses
x = 0.1;
xa = [x, zeros(1,2)];
Pa = diag(ones(1,3));
Pa(1,1) = 0.1;

Wm = [lambda/(L+lambda), 0.5/(L+lambda)+zeros(1,2*L)];      %Mean weights

Wc=Wm;
Wc(1)= Wc(1) + (1-alpha^2+beta);                            %Covariance weights

for t = 2:n-1 

    vol(t-1) = sqrt(x); 

    Xa(:,1) = xa';

    for i = 1:L
        for j = 1:L
            if i == j
                if Pa(i,j) < 0.000001
                    Pa(i,j) = 0.0000000001;
                end
            else
                if Pa(i,j) < 0.000001
                    Pa(i,j) = 0;
                end
            end
        end
    end

    [~,eigPa] = chol(Pa);
    while eigPa ~= 0    %Checking for positive definiteness
        Pa = Pa + 0.001*diag(ones(1,length(Pa)));
        [~,eigPa] = chol(Pa);
    end

    Pa_chol = chol(Pa);

    %% SIGMA POINTS

    for l = 2:L+1
        for i = 1:L
            Xa(i,l) = xa(i) + (sqrt(L+lambda)*Pa_chol(i,l-1));
        end
    end

    for l = 2+L:2*L+1
        for i = 1:L
            Xa(i,l) = xa(i) - (sqrt(L+lambda)*Pa_chol(i,l-L-1));
        end
    end

    % Cholesky
    for l = 1:2*L+1
        if Xa(1,l) < 0
            Xa(1,l) = 0.0001; %Ensuring positivity
        end

        X(l) = Xa(1,l) + (kappa*theta - (mu*rho*sigma) - (kappa - 0.5*rho*sigma)*Xa(1,l))*dt...
            + rho*sigma*(Stock(t) - Stock(t-1))...
            + sigma*sqrt((1-(rho^2))*dt*Xa(1,l))*Xa(2,l);
    end

    x1 = 0;
    for l = 1:2*L+1
        x1 = x1 + Wm(l)*X(l);
    end

    P1 = 0;
    for l = 1:2*L+1
        P1 = P1 + Wc(l)*((X(l)-x1)^2);
    end

    yhat = 0;
    for l = 1:2*L+1
        if X(l) < 0
            X(l) = 0.00001;
        end
        Y(l) = Stock(t) + (mu - 0.5*X(l))*dt + sqrt(X(l))*Xa(3,l);   

        yhat = yhat + Wm(l)*Y(l);
    end

   %% MEASUREMENT UPDATE

    %Transform cross-covariance
    Pyy = 0;
    for l = 1:2*L+1
        Pyy = Pyy + Wc(l)*((Y(l)-yhat)^2);
    end

    %Transform covariance
    Pxy = 0;
    for l = 1:2*L+1
        Pxy = Pxy + Wc(l)*(X(l)-x1)*(Y(l)-yhat);
    end

    K = Pxy/Pyy;                %Kalman Gain

    u(t) = Stock(t+1) - yhat;    %u(count)   %Mean observation errors
    v(t) = Pyy;                 %v(count)   %Variance of observation errors
    estimates(t+1) = yhat;      %estimates(count)

    x = x1 + K*(Stock(t+1)-yhat);    %State update
    P = P1 - ((K^2)*Pyy);           %Covariance update

    xa(1) = x;
    Pa(1,1) = P;

    if x < 0
        x=0.0001;
    end

    Pa(2,1) = 0;
    Pa(1,2) = 0;
    Pa(2,3) = 0;
    Pa(3,2) = 0;
    Pa(1,3) = 0;
    Pa(3,1) = 0;

end

scaled_vol = rescale(vol,min(data.Volatility),max(data.Volatility));
writematrix(transpose(scaled_vol),"UKF_Vol_3.csv");
smoothed_vol = smoothdata(scaled_vol,"gaussian",5);
writematrix(transpose(smoothed_vol),"UKF_Vol_4.csv");

%% Run only to see the graph (not used inside R code)
%figure(1)
%plot((1:1:length(vol)+1),[data.Volatility(1),smoothed_vol]); drawnow;
%title('UKF','fontSize',15,'fontWeight','bold')
%legend('Filtered vol')
