% nonlinear plate model
% Riks method

% clear all;

%load data from previous run

% material parameters
E2D = 349;      % 2D modulus in N/m
v = 0.149;
D = 0.240e-18;  % bending modulus in N*m
%D_G = 0.112e-18;    
%phi = 75*(1-v^2)/8/(23+18*v-3*v^2);
delta0 = 0.6;   % unit: nm
gamma = 0.7;    % unit: J/m^2
k1 = 7.5e18;   % unit: N/m^3
%k2 = 0;
h = sqrt(12*(1-v^2)*D/E2D)*1e9;     % effective thickness in nm

% dimensions
a = 1500;   % radius in nm
b = 250;    % island radius in nm

% normalization 
a = a/h;
b = b/h;
beta = 9*gamma*(h*1e-9)^3/D/2/(delta0*1e-9);
delta0 = delta0/h;
k1 = k1*(h*1e-9)^4/D;
%k2 = k2*(h*1e-9)^5/D;

% numerical parameters
n = 3000;       % number of nodes
n_b = round(b/a*n);

if(n_b == n)
    n_b = n-1;
end

dr = a/n;

al = 5e-1;      % arc length in p-h
toleranceb = 1e-9;      % tolerance for residual 
toleranceR = 5e2;       % tolerance for relative correction 
ntotal = 1;             % total number of steps


% specify the direction of the arc length
al_p = al/sqrt(2);
al_h = al/sqrt(2);
al_v = al/sqrt(2);

%Normalization factor for arc length:  0 represents dimensional; 1 represents dimensionless 
scale = 1;           

%Normalization factor for volume, because volume number is too large in nm^3 unit
v_scale = 1e10;

%Reduce correction vector to avoid blowing up
d_redu = 2;

%Choose control pattern
cont_p = 1;

% cont_p = 1    Pressure control, start with specified pressure
% cont_p = 2    Riks method with pressure and height
% cont_p = 3    Riks method with pressure and volume


%Specify pressure for pressure control
if(cont_p == 1)
    
    large = 0;            %large = 1     Large pressure, using nonlinear plate without vdw solution as initial approximation
                          %large = 0     Small pressure, using zero deflection as initial approxiamtion
                   

    p_p = 500;           %Specify pressure
    
    z = zeros(n-1,1);
    w = zeros(n,1);
    u = zeros(n-1,1);
    

% find the last data point from the previous run for Riks method
else
for i=length(p):-1:1
    if(p(i)~=0)
        break;
    end
end
p_p = p(i);

if(cont_p == 2)
hc_p = height(i);
end

if(cont_p == 3)
volume_p = volume(i);
end

%If starting from origin, first creating variable p, height or volunme and setting it 0
if(p_p == 0)    
    z = zeros(n-1,1);
    w = zeros(n,1);
    u = zeros(n-1,1);
end

end


fzl = zeros(n-1,n-1);
fzn = zeros(n-1,n-1);
fv = zeros(n-1,n-1);
vw = zeros(n-1,n-1);
wz = zeros(n-1,n-1);
% fz = zeros(n-1,n-1);
fu = zeros(n-1,n-1);
fp = zeros(n-1,1);
gz = zeros(n-1,n-1);
gu = zeros(n-1,n-1);
phihz = zeros(1,n-1);
phivw = zeros(1,n-1);

f = zeros(n-1,1);
g = zeros(n-1,1);
% delta = zeros(2*n,1);
% e = zeros(2*n,2*n);
bb = ones(2*n,1);
R = ones(2*n,1);
m = zeros(2*n,1);

wf = zeros(n+1,1);

vdw = zeros(n-1,1);
vdwforce = zeros(n,1);

height = zeros(100,1);
p = zeros(100,1);
volume = zeros(100,1);
result = zeros(100,3);





%winitial = zeros(n+1,1);
%epsilon_r = zeros(n-1,1);
%epsilon_z = zeros(n-1,1);
%z_g = zeros(n-1,1);
%Uall = zeros(300,1);


%Compose Jaccobi matrix (Constant)
   %df/dz linear part
      fzl(1,1) = -1/a^2*(2*n^2+n^2);
      fzl(1,2) = 1/a^2*(n^2+n^2/2);
      
      for k = 2:(n-2)     
      fzl(k,k-1) = 1/a^2*(n^2-n^2/2/k);
      fzl(k,k) = -1/a^2*(2*n^2+n^2/k^2);
      fzl(k,k+1) = 1/a^2*(n^2+n^2/2/k);
      end
  
      fzl(n-1,n-2) = 1/a^2*(n^2-n^2/2/(n-1));
      fzl(n-1,n-1) = -1/a^2*(2*n^2+n^2/(n-1)^2);

    %df/dz vdw part
      %df/dvdw (linear)
      for i = 1:n-1
          for j = 1:i-1
              fv(i,j) = a/n/i*j;
          end
        fv(i,i) = a/n*1/2;
      end

    %dw/dz (linear)
      for i = 1:n-1
          wz(i,i) = -1/2*a/n;
          for j = i+1:n-1
              wz(i,j) = -a/n;
          end 
      end

    %df/dp part (linear)
      for k = 1:(n-1)     
      fp(k,1) = - 1/2*k*a/n;
      end

     %dg/du part (linear)
      gu(1,1) = -1/a^2*(2*n^2+n^2);
      gu(1,2) = 1/a^2*(n^2+n^2/2);
      
      for k = 2:(n-2)       
      gu(k,k-1) = 1/a^2*(n^2-n^2/2/k);
      gu(k,k) = -1/a^2*(2*n^2+n^2/k^2);
      gu(k,k+1) = 1/a^2*(n^2+n^2/2/k);
      end
          
      gu(n-1,n-2) = 1/a^2*(n^2-n^2/2/(n-1));
      gu(n-1,n-1) = -1/a^2*(2*n^2+n^2/(n-1)^2);
      
      
%Pressure control with larger pressure, using nonlinear plate solution without vdw as initial approximation     
if(cont_p == 1 && large == 1)
    
    q = p_p*(h*1e-9)^3/D;
    
    z = 1./16.*q.*a^3.*(-1+1./n.^2.*(1:n-1).^2)./n.*(1:n-1);                  %Using the linear solution to initialize
    z = z';
          
while(norm(bb)>toleranceb||norm(R)>toleranceR)

  %Residue vector
     
     f(1) = 1/a^2*(n^2+n^2/2)*z(2) - 1/a^2*(2*n^2+n^2)*z(1) - 6*n/a*z(1)*u(2) - 12*v*n/a*u(1)*z(1) - 6*z(1)^3 - 1/2*q*a/n;
     g(1) = 1/a^2*(n^2+n^2/2)*u(2) - 1/a^2*(2*n^2+n^2)*u(1) + n/a*(1-v)/2*z(1)^2 + n/2/a*z(1)*z(2);
     
     for k = 2:(n-2)     
      f(k) = 1/a^2*(n^2+n^2/2/k)*z(k+1) - 1/a^2*(2*n^2+n^2/k^2)*z(k) + 1/a^2*(n^2-n^2/2/k)*z(k-1) - 6*n/a*z(k)*(u(k+1)-u(k-1)) - 12*v*n/k/a*u(k)*z(k) - 6*z(k)^3 - 1/2*q*k*a/n;
      g(k) = 1/a^2*(n^2+n^2/2/k)*u(k+1) - 1/a^2*(2*n^2+n^2/k^2)*u(k) + 1/a^2*(n^2-n^2/2/k)*u(k-1) + n/k/a*(1-v)/2*z(k)^2 + n/2/a*z(k)*(z(k+1)-z(k-1));
     end
        
     f(n-1) = - 1/a^2*(2*n^2+n^2/(n-1)^2)*z(n-1) + 1/a^2*(n^2-n^2/2/(n-1))*z(n-2) - 6*n/a*z(n-1)*(-u(n-2)) - 12*v*n/(n-1)/a*u(n-1)*z(n-1) - 6*z(n-1)^3 - 1/2*q*(n-1)*a/n;
     g(n-1) = - 1/a^2*(2*n^2+n^2/(n-1)^2)*u(n-1) + 1/a^2*(n^2-n^2/2/(n-1))*u(n-2) + n/(n-1)/a*(1-v)/2*z(n-1)^2 + n/2/a*z(n-1)*(-z(n-2));
   
 %Compose Jaccobi matrix
    %df/dz nonlinear part
      fzn(1,1) = -6*n/a*u(2) - 12*v*n/a*u(1) - 18*z(1)^2;
      
      for k = 2:(n-2)     
      fzn(k,k) = -6*n/a*(u(k+1)-u(k-1)) - 12*v*n/k/a*u(k) - 18*z(k)^2;
      end
       
      fzn(n-1,n-1) = -6*n/a*(-u(n-2)) - 12*v*n/(n-1)/a*u(n-1) - 18*z(n-1)^2;

    %df/du part (nonlinear)
      fu(1,1) = -12*v*n/a*z(1);
      fu(1,2) = -6*n/a*z(1);
   
      for k = 2:(n-2)     
      fu(k,k-1) = 6*n/a*z(k);
      fu(k,k) = -12*v*n/k/a*z(k);
      fu(k,k+1) = -6*n/a*z(k);
      end
      
      fu(n-1,n-2) = 6*n/a*z(n-1);
      fu(n-1,n-1) = -12*v*n/(n-1)/a*z(n-1);
    
   %dg/dz part (nonlinear)     
      gz(1,1) = n/a*(1-v)*z(1) + n/2/a*z(2);
      gz(1,2) = n/2/a*z(1);
      
      for k = 2:(n-2)     
      gz(k,k-1) = -n/2/a*z(k);
      gz(k,k) = n/k/a*(1-v)*z(k) + n/2/a*(z(k+1)-z(k-1));
      gz(k,k+1) = n/2/a*z(k);
      end
           
      gz(n-1,n-2) = -n/2/a*z(n-1);
      gz(n-1,n-1) = n/(n-1)/a*(1-v)*z(n-1) + n/2/a*(-z(n-2));
      
  
   %Assembly   
    fz = fzl + fzn;  
         
     e = [fz fu;
          gz gu];
      
   m(1:n-1) = z(1:n-1);
   m(n:2*n-2) = u(1:n-1);
   
   while(length(m) ~= 2*(n-1)) 
        m(2*n-1) = [];
   end
   
   bb(1:n-1) = f(1:n-1);
   bb(n:2*n-2) = g(1:n-1);
   
   while(length(bb) ~= 2*(n-1)) 
        bb(2*n-1) = [];
   end
   
   delta = -e\bb;
   
   m = m + delta;
   
   z(1:n-1) = m(1:n-1);
   u(1:n-1) = m(n:2*n-2);
   
end

w(n) = 0;
w(n-1) = 0 - 1/2*(0 + z(n-1))*a/n;

for i = n-2:-1:1
    w(i) = w(i+1) - 1/2*(z(i)+z(i+1))*a/n;
end
   w0 = w(1) - 1/2*z(1)*a/n;  
  
end
      
      

      
  if(cont_p == 2)    
    %dphi_h/dz part (linear)
      for i=1:n-1
       phihz(1,i) = -a/n;
      end
      
    %dphi_h/dh part (linear)
      phihh = -1;
  end
     
   if(cont_p == 3)
     %dphi_v/dz part (linear)
        %dphi_v/dw
          for i=1:n-1
            phivw(1,i) = 2*pi*i*(a/n)^2;
          end
      phivz = phivw*wz/v_scale;    
      
    %dphi_v/dv part (linear)
      phivv = -1/v_scale;
   end
      
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



%Iteration
for mf=1:1:ntotal

    mf

if(cont_p == 1)    
    q = p_p*(h*1e-9)^3/D;
    
else
    if(mf ~= 1)
    p_p = p(mf-1);
    end
    p_p = p_p*(h*1e-9)^3/D;
    q = p_p + al_p*(p_p)^(scale);
        
    
    if(cont_p == 2)
        if(mf ~= 1)
        hc_p = height(mf-1);
        end
      hc_p = hc_p/h;
      hc = hc_p + al_h*(hc_p)^(scale);
    end
    
    if(cont_p == 3)
        if(mf ~= 1)
        volume_p = volume(mf-1);
        end
      volume_p = volume_p/h^3;
      volumec = volume_p + al_v*(volume_p)^(scale);
    end
    
end
  
bb = ones(2*n,1);  
R = ones(2*n,1);


while(norm(bb)>toleranceb||norm(R)>toleranceR)

%     e = zeros(2*n,2*n);
     
     % calculating the vdw force at each node
     for k = 1:n_b
         if(w(k)>=0)
             vdwforce(k) = beta*((delta0/(delta0+w(k)))^4-(delta0/(delta0+w(k)))^10);
         else
             vdwforce(k) = k1*w(k);
         end
     end

     % integrate the vdw force
%     vdw = zeros(n,1);
     vdw(1) = 0.5*a^2/n^2*vdwforce(1);
     for k = 2:n_b
        vdw(k) = vdw(k-1) + 0.5*dr^2*((k-1)*vdwforce(k-1) + k*vdwforce(k));
     end
     
     for k = n_b+1:n-1
        vdw(k) = vdw(n_b);
     end
     
     %integrate theta to obtain deflection
%      w(n) = 0;
%      w(n-1) = 0 - 1/2*(0 + z(n-1))*a/n;
%      for i = n-2:-1:1
%        w(i) = w(i+1) - 1/2*(z(i)+z(i+1))*a/n;
%      end
%      w0 = w(1) - 1/2*z(1)*a/n;
     
     
     %Residue vector
     
     f(1) = 1/a^2*(n^2+n^2/2)*z(2) - 1/a^2*(2*n^2+n^2)*z(1) - 6*n/a*z(1)*u(2) - 12*v*n/a*u(1)*z(1) - 6*z(1)^3 - 1/2*q*a/n + n/a*vdw(1);
     g(1) = 1/a^2*(n^2+n^2/2)*u(2) - 1/a^2*(2*n^2+n^2)*u(1) + n/a*(1-v)/2*z(1)^2 + n/2/a*z(1)*z(2);
     
     for k = 2:(n-2)     
      f(k) = 1/a^2*(n^2+n^2/2/k)*z(k+1) - 1/a^2*(2*n^2+n^2/k^2)*z(k) + 1/a^2*(n^2-n^2/2/k)*z(k-1) - 6*n/a*z(k)*(u(k+1)-u(k-1)) - 12*v*n/k/a*u(k)*z(k) - 6*z(k)^3 - 1/2*q*k*a/n + n/k/a*vdw(k);
      g(k) = 1/a^2*(n^2+n^2/2/k)*u(k+1) - 1/a^2*(2*n^2+n^2/k^2)*u(k) + 1/a^2*(n^2-n^2/2/k)*u(k-1) + n/k/a*(1-v)/2*z(k)^2 + n/2/a*z(k)*(z(k+1)-z(k-1));
     end
        
     f(n-1) = - 1/a^2*(2*n^2+n^2/(n-1)^2)*z(n-1) + 1/a^2*(n^2-n^2/2/(n-1))*z(n-2) - 6*n/a*z(n-1)*(-u(n-2)) - 12*v*n/(n-1)/a*u(n-1)*z(n-1) - 6*z(n-1)^3 - 1/2*q*(n-1)*a/n + n/(n-1)/a*vdw(n-1);
     g(n-1) = - 1/a^2*(2*n^2+n^2/(n-1)^2)*u(n-1) + 1/a^2*(n^2-n^2/2/(n-1))*u(n-2) + n/(n-1)/a*(1-v)/2*z(n-1)^2 + n/2/a*z(n-1)*(-z(n-2));
   
     
     if(cont_p == 2)
       phi_h = w0 - hc;
     
       phi_al = (q - p_p)^2/p_p^(2*scale) + (hc - hc_p)^2/hc_p^(2*scale) - al^2;
     end
     
     if(cont_p == 3)
        phi_v = (v0 - volumec)/v_scale;
        
        phi_al = (q - p_p)^2/p_p^(2*scale) + (volumec - volume_p)^2/volume_p^(2*scale) - al^2;
     end
     
     
%------------------------------------------------------------------------------------------
     
%Compose Jaccobi matrix
    %df/dz nonlinear part
      fzn(1,1) = -6*n/a*u(2) - 12*v*n/a*u(1) - 18*z(1)^2;
      
      for k = 2:(n-2)     
      fzn(k,k) = -6*n/a*(u(k+1)-u(k-1)) - 12*v*n/k/a*u(k) - 18*z(k)^2;
      end
       
      fzn(n-1,n-1) = -6*n/a*(-u(n-2)) - 12*v*n/(n-1)/a*u(n-1) - 18*z(n-1)^2;
 
    %df/dz vdw part
      %dvdw/dw (nonlinear)
      for i = 1:n_b
          if(w(i)>=0)
                vw(i,i) = beta*(10*delta0^10/(delta0 + w(i))^11 - 4*delta0^4/(delta0 + w(i))^5);
          else
                vw(i,i) = k1;
          end
      end
      
      fzv = fv*vw*wz;    
      
    %df/du part (nonlinear)
      fu(1,1) = -12*v*n/a*z(1);
      fu(1,2) = -6*n/a*z(1);
   
      for k = 2:(n-2)     
      fu(k,k-1) = 6*n/a*z(k);
      fu(k,k) = -12*v*n/k/a*z(k);
      fu(k,k+1) = -6*n/a*z(k);
      end
      
      fu(n-1,n-2) = 6*n/a*z(n-1);
      fu(n-1,n-1) = -12*v*n/(n-1)/a*z(n-1);
      

      
    %dg/dz part (nonlinear)     
      gz(1,1) = n/a*(1-v)*z(1) + n/2/a*z(2);
      gz(1,2) = n/2/a*z(1);
      
      for k = 2:(n-2)     
      gz(k,k-1) = -n/2/a*z(k);
      gz(k,k) = n/k/a*(1-v)*z(k) + n/2/a*(z(k+1)-z(k-1));
      gz(k,k+1) = n/2/a*z(k);
      end
           
      gz(n-1,n-2) = -n/2/a*z(n-1);
      gz(n-1,n-1) = n/(n-1)/a*(1-v)*z(n-1) + n/2/a*(-z(n-2));
      

       
    
    if(cont_p == 2)  
    %dphi_al/dp part (nonlinear)  
      phialp = 2*(q - p_p)/p_p^(2*scale);
    
    %dphi_al/dh part (nonlinear)
      phialh = 2*(hc - hc_p)/hc_p^(2*scale);
    end
    
    if(cont_p == 3)  
    %dphi_al/dp part (nonlinear)  
      phialp = 2*(q - p_p)/p_p^(2*scale);
    
    %dphi_al/dv part (nonlinear)
      phialv = 2*(volumec - volume_p)/volume_p^(2*scale);
    end
    
    
      
      

    %Compose zero vector
        z_n1 = zeros(n-1,1);
        z_1n = zeros(1,n-1);
        z_11 = 0;
        
      fz = fzl + fzn + fzv;  
        
    if(cont_p == 1) 
      e = [fz fu;
           gz gu];
    end
      
      
    if(cont_p == 2) 
      e = [fz fu fp z_n1;
           gz gu z_n1 z_n1;
           phihz z_1n z_11 phihh;
           z_1n z_1n phialp phialh];
    end
    
    if(cont_p == 3) 
      e = [fz fu fp z_n1;
           gz gu z_n1 z_n1;
           phivz z_1n z_11 phivv;
           z_1n z_1n phialp phialv];
    end
      
%Residul vector 
   m(1:n-1) = z(1:n-1);
   m(n:2*n-2) = u(1:n-1);
   
   
   if(cont_p == 1) 
      while(length(m) ~= 2*(n-1)) 
        m(2*n-1) = [];
      end
   end
   
   if(cont_p == 2) 
      m(2*n-1) = q;
      m(2*n) = hc;
   end

   if(cont_p == 3) 
      m(2*n-1) = q;
      m(2*n) = volumec;
   end
   
   bb(1:n-1) = f(1:n-1);
   bb(n:2*n-2) = g(1:n-1);
   
   if(cont_p == 1) 
      while(length(bb) ~= 2*(n-1)) 
        bb(2*n-1) = [];
      end
   end
   
   if(cont_p == 2) 
      bb(2*n-1) = phi_h;
      bb(2*n) = phi_al; 
   end
   
   if(cont_p == 3) 
      bb(2*n-1) = phi_v;
      bb(2*n) = phi_al; 
   end  
   
   delta = -e\bb;
   
   R = delta./m;
   fprintf('norm R = %1.5e\n',norm(R))
   
   if(norm(R)>1e0)
       delta = delta/d_redu;
   end
   
   m = m + delta;
   
   z(1:n-1) = m(1:n-1);
   u(1:n-1) = m(n:2*n-2);
   
   
   if(cont_p == 2) 
       q = m(2*n-1);
       hc = m(2*n); 
   end
   
   if(cont_p == 3) 
       q = m(2*n-1);
       volumec = m(2*n); 
   end

   
   w(n) = 0;
   w(n-1) = 0 - 1/2*(0 + z(n-1))*a/n;
   for i = n-2:-1:1
    w(i) = w(i+1) - 1/2*(z(i)+z(i+1))*a/n;
   end
   w0 = w(1) - 1/2*z(1)*a/n;
   
   vtemp = 0;
   for i=1:n-1
   vtemp = vtemp + 2*pi*w(i)*i*a^2/n^2;
   end
   v0 = vtemp;
   
    fprintf('norm bb = %1.5e\n',norm(bb))
    fprintf('\n') 
end 
  
   p(mf) = q*D/(h*1e-9)^3;
   height(mf) = w0*h;
   volume(mf) = v0*h^3;
   
   
x = 0:a*h/n:a*h;
x = x';
for i=n:-1:2;
wf(i) = w(i-1);
end
wf(1) = w0;
wf(n+1) = 0;
ddef = wf*h;

result = [height p volume];

namefig = strcat('profile(p=',num2str(p(mf)),').fig');
% saveas(gcf,'profile_temp.fig')
namemat = strcat('data(p=',num2str(p(mf)),').mat');
save('data_temp.mat','u','v','w','ddef','height','p','volume','result','R','bb','namemat','namefig','x','z','w0','v0','mf');
   
 
% epsilon_r(1) = (u(2)-0)*n/2/a + 1/2*z(1)^2;
% epsilon_z(1) = u(1)*n/a;
% for k =2:n-2
% epsilon_r(k) = (u(k+1)-u(k-1))*n/2/a + 1/2*z(k)^2;
% epsilon_z(k) = u(k)*n/k/a;
% end
% epsilon_r(n-1) = (0-u(n-2))*n/2/a + 1/2*z(n-1)^2;
% epsilon_z(n-1) = u(n-1)*n/(n-1)/a;
%    
% z_g(1) = (z(2)-0)*n/2/a;
% for k =2:n-2
% z_g(k) = (z(k+1)-z(k-1))*n/2/a;
% end
% z_g(n-1) = (0-z(n-2))*n/2/a;  
% r_bend = zeros(n-1,1);
% 
% r_bend(1:n-1,1) = 1:n-1;
% r_bend = r_bend*a/n;   
% 
% u_strech = E2D/2/(1-v^2)*(epsilon_r.^2 + 2*v*epsilon_r.*epsilon_z + epsilon_z.^2);
% 
% u_bend = D/2*(z_g.^2 + z.^2./r_bend.^2 + 2*(D-D_G)/D./r_bend.*z.*z_g)/(h*1e-9)^2;
% 
% u_vdw = -gamma*(3/2*(delta0./(delta0+w)).^3 - 1/2*(delta0./(delta0+w)).^9);
% 
% U_total = 0;
% for i=1:n-1
%     U_total = U_total + 2*pi*(u_strech(i)+u_bend(i)+u_vdw(i))*i*(a*h*1e-9)^2/n^2;
% end
%     U_total = U_total + pi*(u_strech(n-1)+u_bend(n-1)-gamma)*n*(a*h*1e-9)^2/n^2;
%     U_strain = U_total;
%     
%      volume = 0;
% for i=1:n-1
%     volume = volume + 2*pi*(w(i)*h*1e-9)*i*(a*h*1e-9)^2/n^2;
% end
%     volume0(1) = volume;
%     
%     Uall(1) = U_total - q*D/(h*1e-9)^3*volume0(1);
%     Uall(1) = Uall(1)/E2D/(h*1e-9)^2;
%     
%     U_strain = U_strain/E2D/(h*1e-9)^2;
    

end   

x = 0:a*h/n:a*h;
x = x';
for i=n:-1:2;
wf(i) = w(i-1);
end
wf(1) = w0;
wf(n+1) = 0;
ddef = wf*h;

result = [height p volume];

   
% vdwforce = beta.*((delta0./(delta0+w)).^4-(delta0./(delta0+w)).^10);
% 
% vdwforce = vdwforce*D/(h*1e-9)^3*1e-6;

clear delta e epsilon_r epsilon_z f fp fu fv fz fzl fzn fzv g gu gz m phihz r_bend temp_b u_bend u_strech vw winitial wlinear wz z_1n z_g z_n1 zinitial vdw u_vdw pressure
 
namemat = strcat('data(p=',num2str(p(mf)),').mat');
save(namemat);
namefig = strcat('profile(p=',num2str(p(mf)),').fig');
plot(x,ddef)
saveas(gcf,namefig)

% namemat = strcat('data(p=',num2str(p(mf-1)),').mat');
% save(namemat);
% namefig = strcat('profile(p=',num2str(p(mf-1)),').fig');
% plot(x,ddef)
% saveas(gcf,namefig)





