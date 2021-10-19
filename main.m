clear all
close all

%%%%%
% The following script compares different sparse recovery algorithms
% in the case of coherent imaging.
% The methods being compared are
% KM, SQRT_LASSO, StOMP, and TOMP.
% The algorithms try to recover x in 
%   Ax = y
% where A is measurement matrix
%       x is sparse vector
%       y is data vector.

global n_scat

n_scat_range=[8];   % number of scatterers
nsparse=length(n_scat_range);

nreal=1;

noise_level_range=[0.5]; % noise level \delta
nnoise=length(noise_level_range);

phase_new=zeros(nnoise,nsparse);

for jnoise=1:nnoise
    noise_level=noise_level_range(jnoise);
    db_level=floor(-10*log10(noise_level));
    for isp=1:nsparse
        sprintf('isp=%d',isp);
        for nrun=1:nreal
            n_scat=n_scat_range(isp);
                       
            amul=100; % factor that defines array's size in multiple of Lambda
            fmul=2;   % factor the defines bandwidth fraction of 100%
            lmul=200; % factor that defines array's size in multiple of Lambda
            iprint=1; % printing figures with value 1 otherwise use 0
            nm=n_scat;
            
            d = 61;
            nx=d;ny=d; % size of the imaging window
            
            
            n_freq=25; % number of frequencies
            N0=12; Nrec=N0*2+1; % number of receivers on the array
            f_0 =  2400.*10^12/10^5;  % central frequency I am dividing by 10^5 just to avoid numerical zeros
            C0 = 3.*10^8/10^5; % Sound speed
            B = fmul*300.*10^12/10^5; %this is half the bandwidth (corresponding to 100% Bandwidth)
            Lambda = C0/f_0; % the wavelength
            N=n_freq*Nrec;
            noise_size=N;
            nd1=N; % size of data
            K=nx*ny; % size of unknown
            
            model % generates the imaging window
            
            aux = vecnorm(Ac); 
            An = Ac./aux;   % normalize measurement matrix A
            AT = An';
            
            %%%%%
            % The following script applies Kirchhoff migration (KM)
            % The result is x_kir.
            x_kir = An'*y;
            IM_L1=reshape(x_kir,ny,nx)';
            figure(499); PLOT_IMAGE(IM_L1,L);title('Numerics Solution - Kirchhoff '); % the image recovered by KM
            % The end of KM.
            %%%%%
           
           
            tau = thresh(Ac,0,1);    % thresholding parameter \tau
            
            
            %%%%%
            % The following script applies TOMP.
            % The result is x_tomp.
            support_x=[];
            x_tomp = zeros(K,1);
            u = y;
            count = 0;
            tol = 1e-5;
            while count<K
                if norm(u)<tol
                    break;
                else    
                    proxy = abs(Ac(:,1)'*u)/norm(u);
                    count = 0;
                    for i=1:K
                        new_proxy = abs(Ac(:,i)'*u)/norm(u);
                        if new_proxy > tau
                            if new_proxy >= proxy
                                indexmax = i;
                                proxy = new_proxy;
                            end
                        else
                            count = count+1;
                        end
                    end

                    support_x = sort(unique([support_x, indexmax]));
                    u = y - Ac(:,support_x)*((Ac(:,support_x)'*Ac(:,support_x))\(Ac(:,support_x)'*y));
                end    
            end
            x_tomp(support_x) = (Ac(:,support_x)'*Ac(:,support_x))\(Ac(:,support_x)'*y);
            tol = 1e-4;
            x_tomp = max(abs(x_tomp)-tol,0).*sign(x_tomp);  % avoiding machine errors in noiseless case

            IM_L1=reshape(x_tomp,ny,nx)';
            figure(599); PLOT_IMAGE(IM_L1,L);title('Numerics Solution - TOMP'); %  the image recovered by TOMP
            % The end of TOMP.
            %%%%%
            
            
                
            %%%%%
            % The following script applies StOMP using 
            % function stomp(A,y)
            % where A is  measurement matrix, y is data vector.
            % The result is x_stomp.
            [x_stomp,support_x] = stomp(Ac,y);
            IM_L1=reshape(x_stomp,ny,nx)';
            figure(799); PLOT_IMAGE(IM_L1,L);title('Numerics Solution - StOMP'); % the image recovered by StOMP
            % The end of StOMP.
            %%%%%
            
            
            
            %%%%%
            % The following script applies SQRT_LASSO.
            % The result is x_sqrt.
            x_sqrt = zeros(K,1);
            dx = 0.1;
            nsteps = 100;   % number of steps
            eps = 0.001;
            for i=1:nsteps
                update_x = x_sqrt + dx*AT*(y-Ac*x_sqrt)/(norm(y-Ac*x_sqrt)+eps);
                x_sqrt = max(abs(update_x)-tau*dx,0).*sign(update_x);
            end
            IM_L1=reshape(x_sqrt,ny,nx)';
            figure(899); PLOT_IMAGE(IM_L1,L);title('Numerics Solution - SQRTLASSO');    % the image recovered by SQRT_LASSO
            % The end of SQRT_LASSO.
            %%%%%
            
        end
    end
end


%% Vicinity Map
% The following script generates the vicinities 
% S_i = {k: |<a_k,a_i>| > \eta}
% where a_i is a scatterer.

i_scat=sort(i_scat);
strength = 0;
deco = 1/(4*n_scat);    %Vicinity parameter (\eta)

IM = zeros(ix,iy);
IM(IM==0) = NaN;
figure(1111)
hold on
h=pcolor(dom.y-L/Lambda,dom.x,(abs(IM)));
set(h, 'EdgeColor', 'none');
vicinity = ones(K,1);
dx = H_x/Lambda;
dy = H_y/Lambda;
PLOT_IMAGE_back(vicinity,d,dx,dy);  % generating vicinity scene
for i=1:n_scat
    strength = 1;
    vicinity = zeros(K,1);
    vic_strength = [];
    for j=1:K
        if abs(Ac(:,i_scat(i))'*Ac(:,j)) > deco 
            vic_strength = [vic_strength abs(Ac(:,i_scat(i))'*Ac(:,j))];
            vicinity(j) = strength/2;
            if j==i_scat(i)
                vicinity(j) = 2;
            end
        end
    end
   vic_re = reshape(vicinity,ny,nx)';
   PLOT_IMAGE_vic(vicinity,d,dx,dy);    % plotting vicinity points with colors
end
set(gca,'Fontsize',16)
xlabel('range in \lambda_0')
ylabel('cross-range in \lambda_0')
title('Vicinity map');
hold off
%%%%%
%%%%%


%% Plots of solution vectors
% The following script generates the plots of solution vectors
% of different recovery algorithms.

figure(700)
plot(1:K,abs(x_tomp),'k*',1:K,abs(rho),'go')
title('Green circles are the true solution - TOMP')
figure(701)
plot(1:K,abs(x_stomp),'k*',1:K,abs(rho),'go')
title('Green circles are the true solution - StOMP')
figure(702)
plot(1:K,abs(x_sqrt),'k*',1:K,abs(rho),'go')
title('Green circles are the true solution - SQRTLASSO')
figure(703)
plot(1:K,abs(x_kir),'k*',1:K,abs(rho),'go')
title('Green circles are the true solution - Kirchhoff')

%%%%%
%%%%%