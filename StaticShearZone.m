% StaticShearZone.m
% Author: D. Boutelier
% Date: 04/03/2019
% creates synthetic images with a static shear zone 

tic
disp('creating synthetic PIV images with shear zone')
fprintf('\n')

ImSize = [2000 , 2000];     % image size
sizex = ImSize(1);
sizey = ImSize(2);

disp(['image size: ' num2str(sizex) ' x ' num2str(sizey) ' px'])

% preallocate arrays for images
A=zeros(sizey,sizex);

% generate velocity field
x=(-ImSize(1):ImSize(1));

% velocities in pix/frame for both side shear zone
vLeft = 0;      % v left 
vRight = 4;     % v right
disp(['displacement left side: ' num2str(vLeft) ' px'])
disp(['displacement right side: ' num2str(vRight) ' px'])

% position of the limits of shear zone
halfwidthSZ = 256;                      % halfwidth of the shear zone
LeftBound_SZ_X=sizex/2 -halfwidthSZ;    % position of left side of shear zone
RightBound_SZ_X=sizex/2 +halfwidthSZ;   % position of the right side of shear zone
disp(['position left side shear zone: ' num2str(LeftBound_SZ_X) ' px'])
disp(['position right side shear zone: ' num2str(RightBound_SZ_X) ' px'])
fprintf('\n')
fprintf('calculating the velocity field....................')

% width and wavelength
SZ_Width = 2*halfwidthSZ;       % full width of shear zone
w  =4/SZ_Width;                 % wavelength for erf function

% only one velocity component v is calculated - u=0
v=vLeft+0.5*(vRight-vLeft)+0.5*(vRight-vLeft)*erf((w*x)-(w*(LeftBound_SZ_X+SZ_Width/2)));
X=x(x>0);
V=v(x>0);

% make a 2D field
X=repmat(X,ImSize(2),1);
V=repmat(V,ImSize(2),1);

y=[1:ImSize(2)];
Y=repmat(y,ImSize(1),1);
Y=Y';

U=zeros(size(V));
fprintf('done.\n')
fprintf('\n')

%% generate images
% properties of particles
density = 6e-3;         % particle density
Part_D=5;               % particle diameter
Size_Var=4;             % variation diameter px

Sheet_Thick=0.1;        % sheet thickness
noise=0.000;            % noise
Rand_Z=0;               % random Z position percent;

disp(['Particles diameter: ' num2str(Part_D)])

SurfaceArea = sizex*sizey;
N_Particles=round(density*SurfaceArea); disp(['Number of particles: ' num2str(N_Particles)])

% number of images
N_images = 50; 

disp(['Number of images: ' num2str(N_images)]) 

N_i=1;                                          % Image number
randn('state', sum(100*clock));                 %#ok<*RAND>
z0_pre=randn(N_Particles,1);                    % normal distributed sheet intensity
z0=z0_pre*(Rand_Z/200+0.5);

I0=65535*exp(-(Sheet_Thick^2./(0.125*z0.^2)));  % particle intensity
I0(I0>65535)=65535;
I0(I0<0)=0;

randn('state', sum(100*clock));
d=randn(N_Particles,1)/2;                       % particle diameter distribution
d=Part_D+d*Size_Var;
d(d<0)=0;

rand('state', sum(100*clock));
x0=rand(N_Particles,1)*sizex;
y0=rand(N_Particles,1)*sizey;
rd = -8.0 ./ d.^2;
offsety=V;
offsetx=U;


xlimit1=floor(x0-d/2);                          % x min particle extent image1
xlimit2=ceil(x0+d/2);                           % x max particle extent image1
ylimit1=floor(y0-d/2);                          % y min particle extent image1
ylimit2=ceil(y0+d/2);                           % y max particle extent image1
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
                xlimit3(n,1)=floor(x0(n)-d(n)/2-offsetx((y0integer(n)),(x0integer(n)))); 
                xlimit4(n,1)=ceil(x0(n)+d(n)/2-offsetx((y0integer(n)),(x0integer(n)))); 
                ylimit3(n,1)=floor(y0(n)-d(n)/2-offsety((y0integer(n)),(x0integer(n)))); 
                ylimit4(n,1)=ceil(y0(n)+d(n)/2-offsety((y0integer(n)),(x0integer(n)))); 
                
                x1(n)=x0(n)-offsetx((y0integer(n)),(x0integer(n))); 
                y1(n)=y0(n)-offsety((y0integer(n)),(x0integer(n))); 
               
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
                    B(i,j)=B(i,j)+I0(n)*exp((-(j-x0(n)+offsetx(i,j))^2-(i-y0(n)+offsety(i,j))^2)*-rd(n)); 
                end
            end

        end
        fprintf('done.\n')

        B(B>65535)=65535;
        I=imnoise(uint16(B),'gaussian',0,noise);
          
        str = sprintf('%05d',N_i);
        eval(['Frame_' str '=I;']);
        imwrite(eval(['Frame_' str]),['Frame_' str '.tif'],'Compression','none');
        clear I
        clear(['Frame_' str])        
        
    else % image after 2 must be constructed from previous, not from first
        
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
                
            xlimit3(n,1)=floor(x1(n)-d(n)/2-offsetx((y1integer(n)),(x1integer(n)))); 
            xlimit4(n,1)=ceil(x1(n)+d(n)/2-offsetx((y1integer(n)),(x1integer(n)))); 
            ylimit3(n,1)=floor(y1(n)-d(n)/2-offsety((y1integer(n)),(x1integer(n)))); 
            ylimit4(n,1)=ceil(y1(n)+d(n)/2-offsety((y1integer(n)),(x1integer(n)))); 
            
            xnext(n)=x1(n)-offsetx((y1integer(n)),(x1integer(n))); 
            ynext(n)=y1(n)-offsety((y1integer(n)),(x1integer(n))); 

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
                    C(i,j)=C(i,j)+I0(n)*exp((-(j-x1(n)+offsetx(i,j))^2-(i-y1(n)+offsety(i,j))^2)*-rd(n)); 
                end
            end
        end
        fprintf('done.\n')

        C(C>65535)=65535;
        I=imnoise(uint16(C),'gaussian',0,noise);        
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

toc

