////////////////////////////////////////////////////////
// Common functions
////////////////////////////////////////////////////////

function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction


function [Hdot]=HH11(t,x,param)
        Hdot=zeros(5,1);
        Hdot(1)=(1/param(26))*(-param(1)*x(2)*x(3)*(x(1)-param(5)) - param(2)*xinf(x(1),param(10),param(15))*(x(1)-param(6)) - param(3)*x(4)*x(5)*(x(1)-param(6)) - param(4)*(x(1)-param(7)) + I)
        Hdot(2)=(xinf(x(1),param(8),param(13))-x(2))/param(18)
        Hdot(3)=(xinf(x(1),param(9),param(14))-x(3))/param(19)
        Hdot(4)=(xinf(x(1),param(11),param(16))-x(4))/param(20)
        Hdot(5)=(xinf(x(1),param(12),param(17))-x(5))/param(21)
endfunction

function [Hdot]=HH11_Cap_Kir_K_L(t,x,param)
    Hdot=zeros(4,1);
    Hdot(1)=(1/param(22))*(-param(1)*x(2)*(x(1)-param(5)) - param(2)*xinf(x(1),param(9),param(13))*(x(1)-param(6)) - param(3)*x(3)*x(4)*(x(1)-param(6)) - param(4)*(x(1)-param(7)) + I)
    Hdot(2)=(xinf(x(1),param(8),param(12))-x(2))/param(16)
    Hdot(3)=(xinf(x(1),param(10),param(14))-x(3))/param(17)
    Hdot(4)=(xinf(x(1),param(11),param(15))-x(4))/param(18)
endfunction


function [Hdot]=HH22(t,x,param)
    Hdot=zeros(3,1);
    Hdot(1)=(1/param(13))*(-param(1)*x(2)*x(3)*(x(1)-param(3))-param(2)*(x(1)-param(4))+I)
    Hdot(2)=(xinf(x(1),param(5),param(7))-x(2))/param(9)
    Hdot(3)=(xinf(x(1),param(6),param(8))-x(3))/param(10)
endfunction


function [Hdot]=HH22_Kir_K_L(t,x,param)
    Hdot=zeros(3,1);
    Hdot(1)=(1/param(16))*(-param(1)*xinf(x(1),param(6),param(9))*(x(1)-param(4))-param(2)*x(2)*x(3)*(x(1)-param(4))-param(3)*(x(1)-param(5))+I)
    Hdot(2)=(xinf(x(1),param(7),param(10))-x(2))/param(12)
    Hdot(3)=(xinf(x(1),param(8),param(11))-x(3))/param(13)
endfunction

function [Hdot]=HH11_Cat_K_L(t,x,param)
    Hdot=zeros(5,1);
    Hdot(1)=(1/param(23))*(-param(1)*x(2)*x(3)*(x(1)-param(4))-param(2)*x(4)*x(5)*(x(1)-param(5))-param(3)*(x(1)-param(6))+I)
    Hdot(2)=(xinf(x(1),param(7),param(11))-x(2))/param(15)
    Hdot(3)=(xinf(x(1),param(8),param(12))-x(3))/param(16)
    Hdot(4)=(xinf(x(1),param(9),param(13))-x(4))/param(17)
    Hdot(5)=(xinf(x(1),param(10),param(14))-x(5))/param(18)
endfunction


function [Hdot]=HH21(t,x,param)
    Hdot=zeros(4,1);
    Hdot(1)=(1/param(19))*(-param(1)*x(2)*(x(1)-param(4))-param(2)*x(3)*x(4)*(x(1)-param(5))-param(3)*(x(1)-param(6))+I)
    Hdot(2)=(xinf(x(1),param(7),param(10))-x(2))/param(13)
    Hdot(3)=(xinf(x(1),param(8),param(11))-x(3))/param(14)
    Hdot(4)=(xinf(x(1),param(9),param(12))-x(4))/param(15)
endfunction

////////////////////////////////////////////////////////
// Modelling RIM K L
////////////////////////////////////////////////////////

function [Xmin,Xmax]=constRIM_K_L()
    Xmin=[0.1 0.1 -100 -90 -90 -90 1  -30  1  1  0.0001 0.001 0.001];
    Xmax=[50  50  -2   30  -2  -2  30 -1   15 15 0.999  0.999 10];
endfunction

function [model,param]=modelRIM_K_L(param)
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);

    [Xmin,Xmax]=constRIM_K_L();
    D=length(Xmin);
    for i=1:D
        if      param(i)<=Xmin(i) then param(i)=Xmin(i);
        elseif  param(i)> Xmax(i) then param(i)=Xmax(i);
        end
    end
    
    condini = [-38; param(11); param(12)]
    model=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH22); 
        V=x(1,:);
        model(i,:)=V;
    end
endfunction

////////////////////////////////////////////////////////
// Modelling RIM Kir K L
////////////////////////////////////////////////////////

function [Xmin,Xmax]=constRIM_Kir_K_L()
    Xmin=[0.1 0.1 0.1 -100 -90 -90 -90 -90 -30 1  -30  1  1  0.0001 0.001 0.001];
    Xmax=[30  30  30  -2   30  -2  -2  -2  -1  30 -1   15 15 0.999  0.999 10];
endfunction

function [model,param]=modelRIM_Kir_K_L(param)
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    [Xmin,Xmax]=constRIM_Kir_K_L();
    D=length(Xmin);
    for i=1:D
        if      param(i)<=Xmin(i) then param(i)=Xmin(i);
        elseif  param(i)> Xmax(i) then param(i)=Xmax(i);
        end
    end
    
    condini = [-38; param(14); param(15)]
    model=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH22_Kir_K_L); 
        V=x(1,:);
        model(i,:)=V;
    end
endfunction

////////////////////////////////////////////////////////
// Modeling RIM Ca,t K L
////////////////////////////////////////////////////////

function [Xmin,Xmax]=constRIM_Cat_K_L()
    Xmin=[0.1 0.1 0.1 20  -100 -90 -90 -90 -90 -90  1  -30  1  -30 0.0001 0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001 0.001];
    Xmax=[50  50  50  150 -2   30  -2  -2  -2  -2   30 -1   30 -1  15     15     15     15     0.999 0.999 0.999 0.999 10];
endfunction

function [model,param]=modelRIM_Cat_K_L(param)
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    [Xmin,Xmax]=constRIM_Cat_K_L();
    
    D=length(Xmin);
    for i=1:D
        if      param(i)<=Xmin(i) then param(i)=Xmin(i);
        elseif  param(i)> Xmax(i) then param(i)=Xmax(i);
        end
    end
    
    condini = [-38; param(19); param(20); param(21); param(22)]
    model=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH11_Cat_K_L); 
        V=x(1,:);
        model(i,:)=V;
    end
endfunction

////////////////////////////////////////////////////////
// Modelling RIM Ca,p K L
////////////////////////////////////////////////////////

function [Xmin,Xmax]=constRIM_Cap_K_L()
    Xmin=[0.1 0.1 0.1 20  -100 -90 -90 -90 -90  1  1  -30 0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001];
    Xmax=[50  50  50  150 -2   30  -2  -2  -2   30 30 -1  15     15     15     0.999 0.999 0.999 10];
endfunction

function [model,param]=modelRIM_Cap_K_L(param)
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    [Xmin,Xmax]=constRIM_Cap_K_L();
    D=length(Xmin);
    for i=1:D
        if      param(i)<=Xmin(i) then param(i)=Xmin(i);
        elseif  param(i)> Xmax(i) then param(i)=Xmax(i);
        end
    end
    
    condini = [-38; param(16); param(17); param(18)]
    model=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH21); 
        V=x(1,:);
        model(i,:)=V;
    end
endfunction


////////////////////////////////////////////////////////
// Modelling RIM Ca,t Kir K L
////////////////////////////////////////////////////////

function [Xmin,Xmax]=constRIM_Cat_Kir_K_L()
    Xmin=[0.1 0.1 0.1 0.1 20  -100 -90 -90 -90 -90 -90 -90 1  -30 -30 1  -30 0.0001 0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001 0.001];
    Xmax=[50  50  50  50  150 -2   30  -2  -2  -2  -2  -2   30 -1  -1  30 -1  15     15     15     15     0.999 0.999 0.999 0.999 10];
endfunction

function [model,param]=modelRIM_Cat_Kir_K_L(param)
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    [Xmin,Xmax]=constRIM_Cat_Kir_K_L();
    D=length(Xmin);
    for i=1:D
        if      param(i)<=Xmin(i) then param(i)=Xmin(i);
        elseif  param(i)> Xmax(i) then param(i)=Xmax(i);
        end
    end
    
    condini = [-38; param(22); param(23); param(24); param(25)]
    model=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH11); 
        V=x(1,:);
        model(i,:)=V;
    end
endfunction

////////////////////////////////////////////////////////
// Modeling RIM Ca,p Kir K L
////////////////////////////////////////////////////////

function [Xmin,Xmax]=constRIM_Cap_Kir_K_L()
    Xmin=[0.1 0.1 0.1 0.1 20  -100 -90 -90 -90 -90 -90 1  -30 1  -30 0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001];
    Xmax=[50  50  50  50  150 -2   30  -2  -2  -2  -2  30 -1  30 -1  15     15     15     0.999 0.999 0.999 10];
endfunction

function [model,param]=modelRIM_Cap_Kir_K_L(param)
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    [Xmin,Xmax]=constRIM_Cap_Kir_K_L();
    D=length(Xmin);
    for i=1:D
        if      param(i)<=Xmin(i) then param(i)=Xmin(i);
        elseif  param(i)> Xmax(i) then param(i)=Xmax(i);
        end
    end
    
    condini = [-38; param(19); param(20); param(21)]
    model=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH11_Cap_Kir_K_L); 
        V=x(1,:);
        model(i,:)=V;
    end
endfunction

