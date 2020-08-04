syms q h ep H V U t mul mur kB T fl fr z sigma1 sigma2 k a 


%pure system - 3 percent, 1 percent, 0.5 percent
%mixed system- 0.5-3-0.5 percent 

q=1.602e-19;%charge of electron
h=6.636e-34;%planck's constant

m=9.1e-31;%mass of electron

T=300;%room temperature


mul=4.3*q+q*V; %fermi-energy of left with bias V
mur=4.3*q; %fermi-energy of right
 E = sym('E','real'); %energy 

kB=1.38e-23; %boltzmann constant

t_pure=12*0.0005*q; %hopping energy constant t0


%pure system hopping 
%=============================================
t_pure_3_pos=t_pure*exp(10.1257-E/(18.0436*kB*T)); %hopping energy for 3 percent in positive bias
t_pure_3_neg=t_pure*exp(0.4508-E/(232.1742*kB*T));%hopping energy for 3 percent in negative bias

t_pure_1_pos=t_pure*exp(3.0257-E/(61.6523*kB*T));%hopping energy for 1 percent in positive bias
t_pure_1_neg=t_pure*exp(0.358-E/(270.8699*kB*T));%hopping energy for 1 percent in negative bias

t_pure_half_pos=t_pure*exp(2.9257-E/(89*kB*T));%hopping energy for 0.5 percent in positive bias
t_pure_half_neg=t_pure*exp(0.3408-E/(309.5656*kB*T));%hopping energy for 0.5 percent in negative bias

%t_pure=t_pure_3_pos; %active when performing simulation for 3 percent in positive bias
t_pure=t_pure_3_neg; %active when performing simulation for 3 percent in negative bias
%t_pure=t_pure_1_pos;%active when performing simulation for 1 percent in positive bias
%t_pure=t_pure_1_neg; %active when performing simulation for 1 percent in negative bias
%t_pure=t_pure_half_pos;%active when performing simulation for 0.5 percent in positive bias
%t_pure=t_pure_half_neg; %active when performing simulation for 0.5 percent in negative bias


%Mixed system hopping
%====================================
  t_mix_pos=4*0.0185*q; %positive bias hopping constant t0
t_mix_neg=0.0193*8.6*q; %negative bias hopping constant t0
  t_half_neg=0.525*t_mix_neg; %negative bias hopping from 0.5 percent
  t_three_neg= 0.525*t_mix_neg; %negative bias hopping from 3 percent    
%positive bias
   %t_half_pos=0.775*t_mix_pos;
   t_half_pos=0.375*t_mix_pos;
t_three_pos=0;
for v=0.1525:(0.7261-0.1525)/10:0.7261
    %t_three_pos=t_three_pos+exp(-0.2-1.6*v)*t_mix_pos*heaviside(v-V);
    t_three_pos=t_three_pos+exp(-0.7-1.6*v)*t_mix_pos*heaviside(v-V);
end
   %t_three_pos=t_three_pos+0.975*t_mix_pos*heaviside(0.1525-V)+0.25*t_mix_pos;
t_three_pos=t_three_pos+0.975*t_mix_pos*heaviside(0.1525-V)+0.4*t_mix_pos;


%t_three_pos=exp(2.7-5*V)*t_mix_pos+0.275*t_mix_pos;
%t_half_pos=exp(1.7-3*V)*t_mix_pos+0.375*t_mix_pos;

%t_half_mix=t_half_pos; %for positive bias 
%t_three_mix=t_three_pos; %for positive bias 
%t_mix=t_mix_pos;%for positive bias

t_half_mix=t_half_neg; %for negative bias
t_three_mix=t_three_neg; %for negative bias 
t_mix=t_mix_neg;% for negative bias




%rentention energy = concerntration * U
 U_pure_3=0.0179*5*q; % U for 3 percent
 U_pure_1=0.0179*q; % U for 1 percent
 U_pure_half=0.0179*(1/sqrt(2))*q;% U for 0.5 percent
 
 E_g=5.6*q; %energy_gap of PMMA
 
 U_mix=0.0179*2*q;
  U_half_mix=sqrt(1/8)*U_mix;
  U_three_mix=2.5*U_mix;
  %E_g=5.6*q; %commented as already assigned energy gap for PMMA

 
 
 ep_pure=E_g-2*t_pure;% onsite energy for pure system
 ep_mix=E_g-6*t_mix;% onsite energy for mixed system
 
 
 % effective onsite energy = ep_pure+retention energy
 
 a_pure=210e-9; %each length-wise dimension for pure system 
 a_mix=70e-9; %each length-wise dimension for mixed system
 
 

 
 k=abs((2*m*(E-E_g)/(h^2)))^0.5; % momentum from schrodinger equation
%E_max_three=9.398e-19;
%E_max_three_neg=9.385e-19;
%E_max_half=8.977e-19;
%E_max_one=9e-19;


k_three_pos = 4.378e+08;
k_three_neg= 4.1355e+08;
k_one_pos=1.0910e+08;
k_one_neg=1.088e8;
k_half_pos = 4.8960e+07;
k_half_neg=4.8960e+07;


 %Hamiltonian for pure system
 %===================================
 c=3;%concerntration three percent
 %c=1;%concerntration one percent
 %c=0.5;%concerntration half percent
 U_pure=U_pure_3;% activate when concerntration three percent
 %U_pure=U_pure_1;%activate when concerntration one percent
 %U_pure=U_pure_half;%activate when concerntration half percent
 H_pure=ep_pure+c*U_pure;% effective Hamiltonian for pure system = effective onsite energy for pure system
 
 z_pure=exp(-1i*k_three_pos*a_pure);
 sigma1_pure=-t_pure*z_pure; %broadening matrix for pure system at left electrode
 sigma2_pure=-t_pure*z_pure; %broadening matrix for pure system at right electrode
 
 sigma_pure=sigma1_pure+sigma2_pure; %overall broadening matrix for pure system
 
 G_pure=1/(E-H_pure-sigma_pure); %green's function for pure system
 
 tau1_pure=1i*(sigma1_pure-sigma1_pure'); %self-energy matrix at left electrodr for pure system
tau2_pure=1i*(sigma2_pure-sigma2_pure'); %self-energy matrix at right electrode for pure system




%Hamilton for mixed system
%===================================
H_mix=[ep_mix+0.5*U_half_mix -t_half_mix 0;-t_three_mix ep_mix+3*U_three_mix -t_three_mix;0 -t_half_mix ep_mix+0.5*U_half_mix];
 z_mix=exp(-1i*k*a_mix);
 sigma1_mix=[-t_half_mix*z_mix 0 0;0 0 0;0 0 0];
 sigma2_mix=[0 0 0;0 0 0;0 0 -t_half_mix*z_mix];

 sigma_mix=sigma1_mix+sigma2_mix;
 
 G_mix=inv(E*eye(3)-H_mix-sigma_mix); %green's function for mixed system

tau1_mix=1i*(sigma1_mix-sigma1_mix');
 tau2_mix=1i*(sigma2_mix-sigma2_mix');

 
 
 
 
fl=fr+q*(V)/(0.004*kB*T);% difference of fermi distribution of left and right



Trans_pure=trace(tau1_pure*G_pure*tau2_pure*(G_pure'));%transmibility probability for pure system

I_E_pure=q*Trans_pure*(fl-fr)/h; %current as a function of energy I(E)for pure system

I_pure=int(I_E_pure,E,mul,mur); %current I for pure system 



Trans_mix=trace(tau1_mix*G_mix*tau2_mix*(G_mix')); %transmibility probability for mixed system

I_E_mix=q*Trans_mix*(fl-fr)/h;

I_mix=eval(int(I_E_mix,E,mul,mur));



%print
hold on
%pure system
%==============
%fplot(log10(abs(I_pure)),[0 3]) % positive bias for pure system 
fplot(log10(abs(I_pure)),[-3 0]) % negative bias for pure system
%plot(pure_three(:,1),log10(pure_three(:,2))) %actual/experimental result for three percent 
%plot(pure_one(:,1),log10(pure_one(:,2))) %actual/experimental result for three percent
% plot(pure_half(:,1),log10(pure_half(:,2))) %actual/experimental result for three percent
 
 %mixed system
 %===============
%fplot(log10(abs(I_mix)),[0 3]) % positive bias for mixed system 
%fplot(log10(abs(I_mix)),[-3 0]) % negative bias for mixed system
%plot(mix1(:,1),log10(mix1(:,2))) %actual/experimental result for mixed system


%Transmibility probability 
%============================
%fplot(Trans_pure,[mur 2*mur])

