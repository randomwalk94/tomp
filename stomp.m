function [x, support_x, step_count] = stomp(A,dat)

%% Parameters
[N, K] = size(A);

tol = 0.0001;
%% normalizing the columns to 1
aux = vecnorm(A); 
An = A./aux;
AT = An';

%% Thresholding parameter
tau = thresh(A,0,1);

%% TGP

count=0;
support_x=[];
datt=dat;
step_count = 0;
while (count<1)

    x = max(abs(AT*dat) - tau*norm(dat),0);         % Thresholding
    
    nonzero = find(x)'; % find the support at each step
    
    if isempty(nonzero) % if nothing is returned, end the loop
        count=count+1;
    else
        count=0;
    end
    
    support_x = sort(unique([support_x, nonzero]));
    proxy = zeros(K,1);
    proxy(nonzero) = (An(:,nonzero)'*An(:,nonzero))\(An(:,nonzero)'*dat);
    tolx = 1e-6;
    for i=1:K
        if abs(proxy(i))<tolx
            x(i)=0;
        end
    end
    dat = dat - An(:,nonzero)*proxy(nonzero); % project data on the complement of the support
    if norm(dat)<tol % if nothing new is recovered, end the loop
        break;
    end    
%     count = count + 1;
    step_count = step_count + 1;
end
x = zeros(K,1);
x(support_x) = (An(:,support_x)'*An(:,support_x))\(An(:,support_x)'*datt);
tolx = 1e-6;
for i=1:K
    if abs(x(i))<tolx
        x(i)=0;
    end
end


 
end
