%discretization is in steps of resolution
%Lambda *L/a for cross-range
%C0/(2*B) for range
xmul=lmul/amul; ymul=fmul;

% frequencies in the bandwidth
if n_freq==1
    f=f_0;
else
    f=linspace(f_0-B,f_0+B,n_freq);
end
omega = 2*pi*f;

save frequency.mat f_0 Lambda omega C0 n_freq


a = amul*Lambda; %array size in Lambda
L = lmul*Lambda ; %the range in Lambda


y_M(1) = 0;
y_M(2) = L;

Mx = (nx-1)/2;
My = (ny-1)/2;
H_x = Lambda/xmul;
H_y = Lambda/ymul;
ix1 = 1; ix2 = 2*Mx+1;
iy1 = 1; iy2 = 2*My+1;

Oxd =  (y_M(1) + (ix1-(Mx+1))*H_x)/Lambda;
Oyd =  (y_M(2) + (iy1-(My+1))*H_y)/Lambda;

Lx = (ix2-ix1)*H_x/Lambda;
Ly = (iy2-iy1)*H_y/Lambda;

Nldx = Lambda/H_x;
Nldy = Lambda/H_y;

dom = make_domain(Oxd, Oyd, Lx, Ly, Nldx,Nldy);


save dom.mat dom ix1 ix2 iy1 iy2 H_x H_y 


pp= (a/(2*N0))/Lambda ;  % pp : the interelement array distance is pp*Lambda



%the reflectivity
%
rho_scat=(randn(n_scat,1)+sqrt(-1)*randn(n_scat,1))/sqrt(2);
% rho_scat=0.1*randi([1 2],n_scat,1)+sqrt(-1)*0.1*randi([1 2], n_scat,1);
i_scat = randi([1 nx*ny],n_scat,1); 
%
rho=zeros(K,1);
%
for ix = ix1:ix2
    for iy = iy1:iy2
        I_YS = (ix-ix1)*ny + iy-iy1+1;
        xgrid(I_YS) = Lambda*dom.x(ix-ix1+1);
        ygrid(I_YS) = Lambda*dom.y(iy-iy1+1);
    end
end
%



for j=1:n_scat
    rho(i_scat(j))=rho_scat(j);
    y_scat(j,1) = xgrid(i_scat(j)); 
    y_scat(j,2) = ygrid(i_scat(j)); 
end
%

save scat.mat y_scat n_scat y_M C0 Lambda rho_scat i_scat

% source and receivers locations

for k=-N0:N0
    x_s(1,k+N0+1) = k*pp*Lambda ;
    x_s(2,k+N0+1) = 0 ;
end;

Xant = x_s(1,:)/Lambda;
Yant = x_s(2,:)/Lambda;


save array.mat Nrec Xant Yant x_s Lambda C0


%constructing the matrxi A
Ac = zeros(Nrec*n_freq,K);
G = zeros(Nrec,K);
rv = zeros(Nrec,K);


for ix = ix1:ix2
    %  sprintf('ix=%g',ix-ix1+1)
    for iy = iy1:iy2
        I_YS = (ix-ix1)*ny + iy-iy1+1;
        % Calculate the distance from point I_YS to the receivers
        rv(:,I_YS) = sqrt((Lambda*Xant'-Lambda*dom.x(ix-ix1+1)).^2+...
            (Lambda*Yant'-Lambda*dom.y(iy-iy1+1)).^2);
    end
end

I=sqrt(-1);
i=1;
G = exp(-I*omega(i)*rv/C0)/sqrt(N);
Ac = [G];
for i=2:n_freq
    G = exp(-I*omega(i)*rv/C0)/sqrt(N);
    Ac=[Ac; G];
end

%
b = Ac*rho;


% adding noise 
dx=(randn(noise_size,1)+sqrt(-1)*randn(noise_size,1))/sqrt(2);
dx=dx/norm(dx);
nn = noise_level*norm(b)*dx;
y = (b+noise_level*norm(b)*dx); 

%plotting the true solution
IM_true=reshape(rho,ny,nx)';
figure(14); PLOT_IMAGE(IM_true,L);title('True Solution');

%
% The model matrix restricted on the true support 
i=1;
GT = exp(-I*omega(i)*rv(:,i_scat)/C0)/sqrt(N);
AT = [GT];
for i=2:n_freq
    GT = exp(-I*omega(i)*rv(:,i_scat)/C0)/sqrt(N);
    AT=[AT; GT];
end

