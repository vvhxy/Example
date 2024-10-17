clear
clc
%%%%%%%---The known parameters are set as follows
n=4;
M=3;
A_g1=[5.2,0,0,0;
      0,5,0,0;
      0,0,4.9,0;
      0,0,0,5.6];
B_g1=[0.1,-0.2,0.8,-1;
      0.2,-0.4,0.6,0.1;
      0.4,-0.8,0.3,-0.9;
      0.6,0.7,-0.4,1];
C_g1=[1.1,-0.6,-0.8,-0.5;
      1.2,-1.2,1,-1.1;
      1,1,-0.3,-0.9;
      -0.8,1,-0.4,1.1];
d_1=[0 1 1 0];
d_2=[1 0 1 1];
d_3=[0 1 0 1];
D=[d_1;d_2;d_3];
kappa=1.5;
eta_m=0.1;
eta_M=0.3;
L1=-0.6*eye(8);
L2=0.6*eye(8);
rho1=1;
rho2=1;
N=[eye(n),zeros(n)];

%%%%%%%---Define the decision variables
G=sdpvar(n,n); U=sdpvar(n,n); 
q1=sdpvar(n,n); q2=sdpvar(n,n); 
r1=sdpvar(n,n); r2=sdpvar(n,n); 
s1=sdpvar(n,n); s2=sdpvar(n,n); 
w1=sdpvar(n,n); w2=sdpvar(n,n); 
Zg1=sdpvar(1,1); 
R_1=sdpvar(1,1); R_2=sdpvar(1,1); R_3=sdpvar(1,1);
K_1=sdpvar(1,1); K_2=sdpvar(1,1); K_3=sdpvar(1,1);
k1=sdpvar(1,1); k2=sdpvar(1,1); k3=sdpvar(1,1);
O_1=sdpvar(2*n,2*n,'diagonal');
O_2=sdpvar(2*n,2*n,'diagonal');
Omega_1=sdpvar(n,n);
Omega_2=sdpvar(n,n);
Omega_31=sdpvar(n,1); 
Omega_32=sdpvar(n,1);
Omega_33=sdpvar(n,1);
Omega_4=sdpvar(n,n);
Omega_5=sdpvar(n,n);
Omega_6=sdpvar(n,n);
Omega_7=sdpvar(n,n);
Omega_8=sdpvar(n,n);

%%%%%%%---Define the matrices Q_g1,R_g1,S_g1,and W_g1
Q_g1=[q1,zeros(n,n);
     zeros(n,n),q2];
R_g1=[r1,zeros(n,n);
     zeros(n,n),r2];
S_g1=[s1,zeros(n,n);
     zeros(n,n),s2];
W_g1=[w1,zeros(n,n);
     zeros(n,n),w2];

%%%%%%%%%%%%%%---Define the formula (45)
Pi=[G,eye(n);
    eye(n),U];

%%%%%%%---Define the relevant elements in (44)
I_1=[eye(1,1);
     zeros(1,1);
     zeros(1,1)];
I_2=[zeros(1,1);
     eye(1,1);
     zeros(1,1)];
I_3=[zeros(1,1);
     zeros(1,1);
     eye(1,1)];

F_1=2*kappa*R_1-Zg1;
F_2=2*kappa*R_2-Zg1;
F_3=2*kappa*R_3-Zg1;

E_1=[F_1,zeros(1,1); zeros(1,1),F_2];
E_2=[F_1,zeros(1,1); zeros(1,1),F_3];
E_3=[F_2,zeros(1,1); zeros(1,1),F_3];

Xi_1=[A_g1*G+Omega_1,A_g1;Omega_2,U*A_g1];
Xi_2=[zeros(n,n),zeros(n,n);
      Omega_4,Omega_31*I_1'*D+Omega_32*I_2'*D+Omega_33*I_3'*D];
Xi_3=[Omega_5,B_g1;
      Omega_6,U*B_g1];
Xi_4=[Omega_7,C_g1;
      Omega_8,U*C_g1];
Xi_51=[zeros(n,1),zeros(n,1);
       Omega_31,Omega_32];      
Xi_52=[zeros(n,1),zeros(n,1);
       Omega_31,Omega_33];      
Xi_53=[zeros(n,1),zeros(n,1);
       Omega_32,Omega_33];     

H_1=[-Xi_1,Xi_2,zeros(2*n,12*n),Xi_3,Xi_4,Xi_51];

%%%%%%%%%%%%%%---Define an example of the formula (44)
D_11=-Xi_1-Xi_1'+2*kappa*Pi+Q_g1-4*exp(-2*kappa*eta_m)*S_g1-L1*O_1*L2;
D_12=Xi_2;
D_13=-2*exp(-2*kappa*eta_m)*S_g1;
D_16=6*exp(-2*kappa*eta_m)*S_g1;
D_19=Xi_3+0.5*(L1+L2)*O_1;
D_110=Xi_4;
D_111=Xi_51; 
D_14=zeros(2*n,2*n);
D_15=zeros(2*n,2*n);
D_17=zeros(2*n,2*n);
D_18=zeros(2*n,2*n);

D_22=-8*exp(-2*kappa*eta_M)*W_g1;
D_23=-2*exp(-2*kappa*eta_M)*W_g1;
D_24=-2*exp(-2*kappa*eta_M)*W_g1;
D_27=6*exp(-2*kappa*eta_M)*W_g1;
D_28=6*exp(-2*kappa*eta_M)*W_g1;
D_25=zeros(2*n,2*n);
D_26=zeros(2*n,2*n);
D_29=zeros(2*n,2*n);
D_210=zeros(2*n,2*n);
D_211=zeros(2*n,M-1);

D_33=exp(-2*kappa*eta_m)*(R_g1-Q_g1)-4*exp(-2*kappa*eta_m)*S_g1-4*exp(-2*kappa*eta_M)*W_g1;
D_36=6*exp(-2*kappa*eta_m)*S_g1;
D_38=6*exp(-2*kappa*eta_M)*W_g1;
D_34=zeros(2*n,2*n);
D_35=zeros(2*n,2*n);
D_37=zeros(2*n,2*n);
D_39=zeros(2*n,2*n);
D_310=zeros(2*n,2*n);
D_311=zeros(2*n,M-1);

D_44=-exp(-2*kappa*eta_m)*R_g1-4*exp(-2*kappa*eta_M)*W_g1;
D_47=6*exp(-2*kappa*eta_M)*W_g1;
D_45=zeros(2*n,2*n);
D_46=zeros(2*n,2*n);
D_48=zeros(2*n,2*n);
D_49=zeros(2*n,2*n);
D_410=zeros(2*n,2*n);
D_411=zeros(2*n,M-1);

D_55=-L1*O_2*L2;
D_510=0.5*(L1+L2)*O_2;
D_56=zeros(2*n,2*n);
D_57=zeros(2*n,2*n);
D_58=zeros(2*n,2*n);
D_59=zeros(2*n,2*n);
D_511=zeros(2*n,M-1);

D_66=-12*exp(-2*kappa*eta_m)*S_g1;
D_67=zeros(2*n,2*n);
D_68=zeros(2*n,2*n);
D_69=zeros(2*n,2*n);
D_610=zeros(2*n,2*n);
D_611=zeros(2*n,M-1);

D_77=-12*exp(-2*kappa*eta_M)*W_g1;
D_78=zeros(2*n,2*n);
D_79=zeros(2*n,2*n);
D_710=zeros(2*n,2*n);
D_711=zeros(2*n,M-1);

D_88=-12*exp(-2*kappa*eta_M)*W_g1;
D_89=zeros(2*n,2*n);
D_810=zeros(2*n,2*n);
D_811=zeros(2*n,M-1);

D_99=-O_1;
D_910=zeros(2*n,2*n);
D_911=zeros(2*n,M-1);

D_1010=-O_2;
D_1011=zeros(2*n,M-1);

D_1111=E_1;

D1=[D_11,D_12,D_13,D_14,D_15,D_16,D_17,D_18,D_19,D_110,D_111;
   D_12',D_22,D_23,D_24,D_25,D_26,D_27,D_28,D_29,D_210,D_211;
   D_13',D_23',D_33,D_34,D_35,D_36,D_37,D_38,D_39,D_310,D_311;
   D_14',D_24',D_34',D_44,D_45,D_46,D_47,D_48,D_49,D_410,D_411;
   D_15',D_25',D_35',D_45',D_55,D_56,D_57,D_58,D_59,D_510,D_511;
   D_16',D_26',D_36',D_46',D_56',D_66,D_67,D_68,D_69,D_610,D_611;
   D_17',D_27',D_37',D_47',D_57',D_67',D_77,D_78,D_79,D_710,D_711;
   D_18',D_28',D_38',D_48',D_58',D_68',D_78',D_88,D_89,D_810,D_811;
   D_19',D_29',D_39',D_49',D_59',D_69',D_79',D_89',D_99,D_910,D_911;
   D_110',D_210',D_310',D_410',D_510',D_610',D_710',D_810',D_910',D_1010,D_1011;
   D_111',D_211',D_311',D_411',D_511',D_611',D_711',D_811',D_911',D_1011',D_1111];

X1=rho1^2*S_g1-2*rho1*Pi;
X2=rho2^2*W_g1-2*rho2*Pi;
L_1=[eta_m*H_1',(eta_M-eta_m)*H_1',sqrt(eta_M-eta_m)*H_1'*N'*d_1',sqrt(eta_M-eta_m)*H_1'*N'*d_2',sqrt(eta_M-eta_m)*H_1'*N'*d_3'];

L4=[X1,zeros(2*n,2*n+M);
    zeros(2*n,2*n),X2,zeros(2*n,M);
    zeros(1,4*n),-k1,zeros(1,M-1);
    zeros(1,4*n+1),-k2,zeros(1,M-2);
    zeros(1,4*n+M-1),-k3];

Y1=[D1,L_1;
   L_1',L4];

%%%%%%%%%%%%%%---Define an example of the formula (46)
ADD1=[k1,eye(1,1);
      eye(1,1),K_1];

%%%%%%%%%%%%%%---Add the constraint conditions
F=[Y1<=0.01*eye(24*n+2*M-1);
   Pi>=0.01*eye(2*n);
   ADD1>=0.01*eye(2);
   G>=0.01*eye(n);U>=0.01*eye(n);
   Q_g1>=0.01*eye(2*n);R_g1>=0.01*eye(2*n);
   S_g1>=0.01*eye(2*n),W_g1>=0.01*eye(2*n);
   Zg1>=0.01*eye(1);
   O_1>=0.01*eye(2*n);O_2>=0.01*eye(2*n);
   R_1>=0.01*eye(1);R_2>=0.01*eye(1);R_3>=0.01*eye(1);
   K_1>=0.01*eye(1);K_2>=0.01*eye(1);K_3>=0.01*eye(1)];

%%%%%%%---Setup the solver
ops=sdpsettings('solver','sdpt3','verbose',1);

%%%%%%%---Find the feasible solution
solvesdp(F,'',ops)

%%%%%%%---extract the value of the decision variables
G=double(G);U=double(U);
Q_g1=double(Q_g1);R_g1=double(R_g1);S_g1=double(S_g1);W_g1=double(W_g1);
Zg1=double(Zg1);
R_1=double(R_1);R_2=double(R_2);R_3=double(R_3);
K_1=double(K_1);K_2=double(K_2);K_3=double(K_3);
k1=double(k1);k2=double(k2);k1=double(k3);
O_1=double(O_1);O_2=double(O_2);
Omega_1=double(Omega_1);Omega_2=double(Omega_2);
Omega_31=double(Omega_31);Omega_32=double(Omega_32);Omega_33=double(Omega_33);
Omega_4=double(Omega_4);Omega_5=double(Omega_5);
Omega_6=double(Omega_6);Omega_7=double(Omega_7);Omega_8=double(Omega_8);