function RhoOut = EntanglementRotationGate(RhoIn,param,angles)
Rho=RhoIn;

nqbits=length(angles)/3;
angles = reshape(angles,[nqbits,3]);
rabi=param.rabi;

angle1=angles(1);
angle2=angles(2);
angle3=angles(3);

I=eye(4);
X=[0 1;1 0];
Y=[0 -1i;1i 0];
Z=[1 0;0 -1];

Uy1=cos(angle1)*I-1i*sin(angle1)*kron(Y,Y);

Rho=Uy1*Rho*Uy1';
param.entanglement_time=abs(real(angle1))/rabi;
if param.error; Rho=LindbladSelect(Rho,param,'entanglement'); end

Uz=cos(angle2)*I-1i*sin(angle2)*kron(Z,Z);

Rho=Uz*Rho*Uz';
param.entanglement_time=abs(real(angle2))/rabi;
if param.error; Rho=LindbladSelect(Rho,param,'entanglement'); end

Uy2=cos(angle3)*I-1i*sin(angle3)*kron(Y,Y);

Rho=Uy2*Rho*Uy2';
param.entanglement_time=abs(real(angle3))/rabi;
if param.error; Rho=LindbladSelect(Rho,param,'entanglement'); end

RhoOut=Rho;