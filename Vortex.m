% Vortex.m
% Author: D. Boutelier
% Date: 04/03/2019
% creates synthetic images with a vortex 

L = 2000;           % image size
x = (0:L);          % vector x 
y = (0:L);          % vector y

Xc = round(L/2);    % coordinate of vortex center
Yc = round(L/2);

[X,Y] = meshgrid(x,y);
[T,R] = cart2pol(X-Xc,Y-Yc);
[n,m] = size(R);

D = zeros(n,m);
V = zeros(n,m);

Rad = 512+256;      % radius of teh vortex
Vmax = 4;           % magnitude of the velocity at the edge of the vortex
Vm=ones(size(R)).*Vmax;
Vm(R<=Rad) = (Vmax/Rad) .*R(R<=Rad);
Vm(R>Rad) = Vmax./(1+R(R>Rad)-Rad);

U = -Vm.*sin(T);
V = Vm.*cos(T);

% downsize X,Y,U,V for plot
k=100;
XS = X(1:k:end, 1:k:end);
YS = Y(1:k:end, 1:k:end);

US = U(1:k:end, 1:k:end);
VS = V(1:k:end, 1:k:end);

% plot quiver map 
figure(1)
imagesc(X(1,:),Y(:,1),Vm)
hold on
quiver(XS,YS,US,VS,'color', 'red')
daspect([1 1 1])
axis([0 L 0 L])
colorbar()

% plot velocity profile
figure(2)
plot(X(1001,:), V(1001,:), '-')


%% generate images
% properties of particles
density = 6e-3;     % particle density
Part_D = 5;         % particle diameter
Size_Var = 4;       % variation in diameter
Sheet_Thick = 0.1;  % sheet size
noise = 0.000;
Rand_Z = 0;         % percent;

disp(['Particles diameter: ' num2str(Part_D)])

sizex = L;
sizey = L;
A=zeros(sizey,sizex);
SurfaceArea = sizex*sizey;
N_Particles=round(density*SurfaceArea); disp(['Number of particles: ' num2str(N_Particles)])

% number of images
N_images = 10; 

disp(['Number of images: ' num2str(N_images)]) 

N_i=1; %Image number
randn('state', sum(100*clock)); %#ok<*RAND>
z0_pre=randn(N_Particles,1); %normal distributed sheet intensity
z0=z0_pre*(Rand_Z/200+0.5);

I0=65535*exp(-(Sheet_Thick^2./(0.125*z0.^2))); %particle intensity
I0(I0>65535)=65535;
I0(I0<0)=0;

randn('state', sum(100*clock));
d=randn(N_Particles,1)/2; %particle diameter distribution
d=Part_D+d*Size_Var;
d(d<0)=0;

rand('state', sum(100*clock));
x0=rand(N_Particles,1)*sizex;
y0=rand(N_Particles,1)*sizey;
rd = -8.0 ./ d.^2;
offsety=V;
offsetx=U;


xlimit1=floor(x0-d/2); %x min particle extent image1
xlimit2=ceil(x0+d/2); %x max particle extent image1
ylimit1=floor(y0-d/2); %y min particle extent image1
ylimit2=ceil(y0+d/2); %y max particle extent image1
xlimit2(xlimit2>sizex)=sizex;
xlimit1(xlimit1<1)=1;
ylimit2(ylimit2>sizey)=sizey;
ylimit1(ylimit1<1)=1;

%% create image 1

fprintf(['creating image ' num2str(N_i)])
for n=1:N_Particles
    p=n/N_Particles*100;
    pp=mod(p,5);
    if pp == 0
        fprintf('.')
    end
    
    r = rd(n);
    for j=xlimit1(n):xlimit2(n)
        rj = (j-x0(n))^2;
        for i=ylimit1(n):ylimit2(n)
            A(i,j)=A(i,j)+I0(n)*exp((rj+(i-y0(n))^2)*r);
        end
    end
end
 fprintf('done.\n')

A(A>65535)=65535;
I = imnoise(uint16(A),'gaussian',0,noise);
str = sprintf('%05d',N_i);
eval(['Frame_' str '=I;']);
imwrite(eval(['Frame_' str]),['Frame_' str '.tif'],'Compression','none');
clear I
clear(['Frame_' str])

%%
for N_i=2:N_images
    fprintf('\n')

    if N_i==2
        x0integer=round(x0);
        x0integer(x0integer>sizex)=sizex;
        x0integer(x0integer<1)=1;
        y0integer=round(y0);
        y0integer(y0integer>sizey)=sizey;
        y0integer(y0integer<1)=1;

        xlimit3=zeros(N_Particles,1);
        xlimit4=xlimit3;
        ylimit3=xlimit3;
        ylimit4=xlimit3;

       fprintf(['calculating shifts for image ' num2str(N_i)])
       
            for n=1:N_Particles
                p=n/N_Particles*100;
                pp=mod(p,5);
                if pp == 0
                    fprintf('.')
                end
                xlimit3(n,1)=floor(x0(n)-d(n)/2-offsetx((y0integer(n)),(x0integer(n)))); %x min particle extent image2
                xlimit4(n,1)=ceil(x0(n)+d(n)/2-offsetx((y0integer(n)),(x0integer(n)))); %x max particle extent image2
                ylimit3(n,1)=floor(y0(n)-d(n)/2-offsety((y0integer(n)),(x0integer(n)))); %y min particle extent image2
                ylimit4(n,1)=ceil(y0(n)+d(n)/2-offsety((y0integer(n)),(x0integer(n)))); %y max particle extent image2
                
                x1(n)=x0(n)-offsetx((y0integer(n)),(x0integer(n))); %new position of center
                y1(n)=y0(n)-offsety((y0integer(n)),(x0integer(n))); %new position of center
               
            end
            fprintf('done.\n')
        
        x1=x1';
        y1=y1';
        
        xlimit3(xlimit3<1)=1;
        xlimit4(xlimit4>sizex)=sizex;
        ylimit3(ylimit3<1)=1;
        ylimit4(ylimit4>sizey)=sizey;
        
        B=zeros(sizey,sizex);

      fprintf(['creating image ' num2str(N_i)])
        for n=1:N_Particles
            p=n/N_Particles*100;
                pp=mod(p,5);
                if pp == 0
                    fprintf('.')
                end
                
            r = rd(n);
            for j=xlimit3(n):xlimit4(n)
                for i=ylimit3(n):ylimit4(n)
                    B(i,j)=B(i,j)+I0(n)*exp((-(j-x0(n)+offsetx(i,j))^2-(i-y0(n)+offsety(i,j))^2)*-rd(n)); %place particle with gaussian intensity profile
                end
            end

        end
        fprintf('done.\n')

        B(B>65535)=65535;
        I = imnoise(uint16(B),'gaussian',0,noise);  
        str = sprintf('%05d',N_i);
        eval(['Frame_' str '=I;']);
        imwrite(eval(['Frame_' str]),['Frame_' str '.tif'],'Compression','none');
        clear I
        clear(['Frame_' str])       
        
    else
        % image after 2 must be constructed from previous, not from first     
        C=zeros(sizey,sizex);
        
        if N_i >3
            x1=xnext;
            y1=ynext;
        end
        
        x1integer=round(0.5*(xlimit3+xlimit4));
        x1integer(x1integer>sizex)=sizex;
        x1integer(x1integer<1)=1;
        
        y1integer=round(0.5*(ylimit3+ylimit4));
        y1integer(y1integer>sizey)=sizey;
        y1integer(y1integer<1)=1;
        
        fprintf(['calculating shifts for image ' num2str(N_i)])
        
        for n=1:N_Particles
            p=n/N_Particles*100;
                pp=mod(p,5);
                if pp == 0
                    fprintf('.')
                end
                
            xlimit3(n,1)=floor(x1(n)-d(n)/2-offsetx((y1integer(n)),(x1integer(n)))); %x min particle extent image2
            xlimit4(n,1)=ceil(x1(n)+d(n)/2-offsetx((y1integer(n)),(x1integer(n)))); %x max particle extent image2
            ylimit3(n,1)=floor(y1(n)-d(n)/2-offsety((y1integer(n)),(x1integer(n)))); %y min particle extent image2
            ylimit4(n,1)=ceil(y1(n)+d(n)/2-offsety((y1integer(n)),(x1integer(n)))); %y max particle extent image2
            
            xnext(n)=x1(n)-offsetx((y1integer(n)),(x1integer(n))); %new position of center
            ynext(n)=y1(n)-offsety((y1integer(n)),(x1integer(n))); %new position of center
            
        end
        fprintf('done.\n')

        xlimit3(xlimit3<1)=1;
        xlimit4(xlimit4>sizex)=sizex;
        ylimit3(ylimit3<1)=1;
        ylimit4(ylimit4>sizey)=sizey;
        
        
        fprintf(['creating image ' num2str(N_i)])
        for n=1:N_Particles
            p=n/N_Particles*100;
                pp=mod(p,5);
                if pp == 0
                    fprintf('.')
                end
             
            r = rd(n);
            for j=xlimit3(n):xlimit4(n)
                for i=ylimit3(n):ylimit4(n)
                    C(i,j)=C(i,j)+I0(n)*exp((-(j-x1(n)+offsetx(i,j))^2-(i-y1(n)+offsety(i,j))^2)*-rd(n)); % place particle with gaussian intensity profile
                end
            end     
        end
        fprintf('done.\n')

        C(C>65535)=65535;
        I = imnoise(uint16(C),'gaussian',0,noise);        
        str = sprintf('%05d',N_i);
        eval(['Frame_' str '=I;']);
        imwrite(eval(['Frame_' str]),['Frame_' str '.tif'],'Compression','none');
        clear I
        clear(['Frame_' str])
        
        clear C I
        clearvars -regexp ^I_\d{1}$
        clearvars -regexp ^I_\d{2}$
        clearvars -regexp ^I_\d{3}$
        
    end 
end
