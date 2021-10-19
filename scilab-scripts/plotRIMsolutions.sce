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
// Benchmark model RIM K L
////////////////////////////////////////////////////////

function [benchmark]=benchmarkRIM_K_L()
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    param= [13.094735710091133 0.36646134715010581 -99.947221058141579 -19.767684223786659 -88.309771061154748 -89.500840984810097 30 -11.09456936074044 14.515278883157309 1.0784423742897404 0.97016279285488238 0.0051348979242179325 0.0083921416327793902];  
    condini = [-38; param(11); param(12)]
    benchmark=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH22); 
        V=x(1,:);
        benchmark(i,:)=V;
        //plot(t,V,'k')
    end
    //p=get("hdl"); //get handle on current entity (here the polyline entity)
endfunction

function [benchmark]=bestFoundRIM_K_L_is_Cap_K_L()
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    param= [0.10081129099719248 11.663416090251854 0.26542537751326151 20.875869879311576 -97.733183085078736 -35.380208129884657 -22.057780108180836 -75.858942339037583 -84.291424914336531 11.139871872905655 17.259349126500496 -7.9561737110379411 1.0530017463452479 14.563328351015672 1.0170705215369544 0.94754193300821232 0.89777160195738892 0.0066028365187863055 0.0084733124700737822];    
    condini = [-38; param(16); param(17); param(18)]
    benchmark=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH21); 
        V=x(1,:);
        benchmark(i,:)=V;
        //plot(t,V);
    end   
endfunction


function plotRIM_K_L(plotSolution)
    if ~exists("plotSolution","local") then
         plotSolution = %F;
    end
  
    stimuli=[-15:5:35];
    t=linspace(0,50,10000)
    
    benchmark = benchmarkRIM_K_L();
    model     = bestFoundRIM_K_L_is_Cap_K_L();
    
    clf();
    p=gca();
    
    p.foreground=1;
    p.thickness=4;
    p.mark_style=9;
    plot(t,benchmark,'c-')
    f=gcf();
    f.figure_size = [600 730];
    f.children.margins(2)=0.02;

    if plotSolution then
        p.foreground=0;
        //p.thickness=1;
        p.mark_style=9;
        plot(t,model,'k-');
        f=gcf();
        f.children.children(1).children.line_style=3;
        f.children.children(1).children.thickness=2;
        title('$Benchmark: I_K+I_L \quad vs. \quad\\ Model: I_{Ca,p}+I_K+I_L$','fontsize',5)
        //legend(['$Benchmark$';'$Model$'],1)
    else
        title('$Benchmark: \quad I_K+I_L$','fontsize',5)
    end   
        xlabel('$t(10^{-1}s)$','fontsize',5);
        ylabel('$V(t)$','fontsize',5);
endfunction






function [benchmark]=test2RIM_K_L()
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    param= [ 0.4941898 0.2643662 -100 -42.455555 -90 -62.409001 6.6587485 -6.1808329 7.7592828 7.2397041 0.3281542 0.001 0.1539549];  
    condini = [-38; param(11); param(12)]
    benchmark=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH22); 
        V=x(1,:);
        benchmark(i,:)=V;
        plot(t,V);
    end
    xtitle ( "IK+IL" , "t(10^{-1}s)" , "V(t)" );
    //plot(t,benchmark)
endfunction 

////////////////////////////////////////////////////////
// Benchmark model RIM Kir K L
////////////////////////////////////////////////////////
function [benchmark]=benchmarkRIM_Kir_K_L()
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    param= [2.8679946715139488 29.797719589982421 0.41698996862343291 -99.994044426903159 -10.743984489831943 -89.867879606107749 -89.194570055686086 -89.61731180323676 -19.307101474638543 28.387450665777671 -1.9473308473159785 14.951137530423145 1.0046323329491569 0.65697896468124639 0.001 0.04301376492040998];    
    condini = [-38; param(14); param(15)]
    benchmark=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH22_Kir_K_L); 
        V=x(1,:);
        benchmark(i,:)=V;
        //plot(t,V)
    end
    //plot(t,benchmark

endfunction

function [benchmark]=bestFoundRIM_Kir_K_L_is_Kir_K_L()
     stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    param= [2.8589173052423957 29.368810088189136 0.41639554636034126 -99.92184324021737 -10.838110278000391 -89.583378585367399 -90 -89.532956382148598 -19.196183779722077 29.779767179173795 -1.9252760809973912 13.946556719694435 1 0.60834387805019263 0.0010228998482015852 0.04304406049447413];    
    condini = [-38; param(14); param(15)]
    benchmark=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH22_Kir_K_L); 
        V=x(1,:);
        benchmark(i,:)=V;
        //plot(t,V)
    end
    //plot(t,benchmark
    
endfunction

function plotRIM_Kir_K_L(plotSolution)
    if ~exists("plotSolution","local") then
         plotSolution = %F;
    end
  
    stimuli=[-15:5:35];
    t=linspace(0,50,10000)
    
    benchmark = benchmarkRIM_Kir_K_L();
    model     = bestFoundRIM_Kir_K_L_is_Kir_K_L();
    
    clf();
    p=gca();
    
    p.foreground=1;
    p.thickness=4;
    p.mark_style=9;
    plot(t,benchmark,'c-')
    f=gcf();
    f.figure_size = [600 730];
    f.children.margins(2)=0.02;

    if plotSolution then
        p.foreground=0;
        //p.thickness=1;
        p.mark_style=9;
        plot(t,model,'k-');
        f=gcf();
        f.children.children(1).children.line_style=3;
        f.children.children(1).children.thickness=2;
        title('$Benchmark: I_{Kir}+I_K+I_L \quad vs. \quad\\ Model: I_{Kir}+I_K+I_L$','fontsize',5)
        //legend(['$Benchmark$';'$Model$'],1)
    else
        title('$Benchmark: \quad  I_{Kir}+I_K+I_L$','fontsize',5)
    end   
        xlabel('$t(10^{-1}s)$','fontsize',5);
        ylabel('$V(t)$','fontsize',5);
endfunction


////////////////////////////////////////////////////////
// Benchmark model RIM Ca,t K L
////////////////////////////////////////////////////////

function [benchmark]=benchmarkRIM_Cat_K_L()
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    param= [0.53669904914909083 0.92299234478342829 0.28414296151074003 149.83572189405282 -28.045190345588459 -64.665698155987656 -2.7156235767589396 -2.1881402169776791 -21.335689228063984 -2 14.298533238465726 -29.907480746249988 10.753482651510865 -12.820982699435923 0.23965708929474799 0.21184704494511333 0.0001 0.4772492347391229 0.081019189379593903 0.999 0.97212838956993897 0.17986446602982653 0.025263343150056742];    
    condini = [-38; param(19); param(20); param(21); param(22)]
    benchmark=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH11_Cat_K_L); 
        V=x(1,:);
        benchmark(i,:)=V;
        //plot(t,V)
    end
    //plot(t,benchmark)
endfunction

////////////////////////////////////////////////////////
// Benchmark model RIM Ca,p K L
////////////////////////////////////////////////////////

function [benchmark]=benchmarkRIM_Cap_K_L()
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    param= [0.30137168331410008 7.8122076927080659 0.19579161951644042 42.093540294272231 -99.999992085264552 -66.904605592271793 -2.0003055931144305 -52.105337571615209 -90 29.999974543300013 22.651122690805924 -2.0584133051349447 0.21325407674762803 2.3764628162983965 15 0.25814963215638864 0.999 0.0016695469672287615 0.033230807198809789];    
    condini = [-38; param(16); param(17); param(18)]
    benchmark=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH21); 
        V=x(1,:);
        benchmark(i,:)=V;
        //plot(t,V)
    end
    //plot(t,benchmark)
endfunction

function [benchmark]=bestFoundRIM_Cap_K_L_is_Cap_K_L()
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    param= [0.29538475249939389 6.8114415303100859 0.19893062589878874 43.484510532253267 -99.909642681038449 -67.019170956016168 -3.4346174692778852 -54.379751129048806 -89.053126650623199 29.441379229735755 22.655334008987779 -1.8469563919744758 0.20430872662732819 2.4340822704369738 14.396894123918731 0.27863006785178757 0.99724602712054577 0.0019447319352347292 0.032849105268709017];    
    condini = [-38; param(16); param(17); param(18)]
    benchmark=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH21); 
        V=x(1,:);
        benchmark(i,:)=V;
        //plot(t,V)
    end
    //plot(t,benchmark)
endfunction


function plotRIM_Cap_K_L(plotSolution)
    if ~exists("plotSolution","local") then
         plotSolution = %F;
    end
  
    stimuli=[-15:5:35];
    t=linspace(0,50,10000)
    
    benchmark = benchmarkRIM_Cap_K_L();
    model     = bestFoundRIM_Cap_K_L_is_Cap_K_L();
    
    clf();
    p=gca();
    
    p.foreground=1;
    p.thickness=4;
    p.mark_style=9;
    plot(t,benchmark,'c-');
    f=gcf();
    f.figure_size = [600 730];
    f.children.margins(2)=0.02;

    if plotSolution then
        p.foreground=0;
        //p.thickness=1;
        p.mark_style=9;
        plot(t,model,'k-');
        f=gcf();
        f.children.children(1).children.line_style=3;
        f.children.children(1).children.thickness=2;
        title('$Benchmark: I_{Ca,p}+I_K+I_L \quad vs. \quad\\ Model:  I_{Ca,p}+I_K+I_L$','fontsize',5);
        //legend(['$Benchmark$';'$Model$'],1)
    else
        title('$Benchmark: \quad I_{Ca,p}+I_K+I_L$','fontsize',5)
    end   
        xlabel('$t(10^{-1}s)$','fontsize',5);
        ylabel('$V(t)$','fontsize',5);
endfunction


////////////////////////////////////////////////////////
// Benchmark model RIM Ca,t Kir K L
////////////////////////////////////////////////////////

function [benchmark]=benchmarkRIM_Cat_Kir_K_L()
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    param= [0.38550371762530544 0.29531094063772839 0.19854213426212194 0.2578929741463718 149.45661754386518 -99.969585899723697 -70.537503916206916 -2 -2.0567413807496124 -89.994240155431228 -19.181899518139403 -23.545959947461395 22.070594527941161 -29.990934263086125 -1 7.9271723244942036 -12.466427293206401 0.22177919739425361 14.965758431051139 0.243333358691743 14.624032628968939 0.18498125974159876 0.76939308081612934 0.30346851721117168 0.75600198625208725 0.025955842785910805];    
    condini = [-38; param(22); param(23); param(24); param(25)]
    benchmark=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH11); 
        V=x(1,:);
        benchmark(i,:)=V;
        //plot(t,V)
    end
    //plot(t,benchmark)
endfunction

////////////////////////////////////////////////////////
// Benchmark model RIM Ca,p Kir K L
////////////////////////////////////////////////////////

function [benchmark]=benchmarkRIM_Cap_Kir_K_L()
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    param= [0.2401245091237183 0.33243551676834415 0.12728412756569329 0.2850521829435072 105.28653202706991 -100 -81.330284879273279 -21.039865777213976 -89.990094699261817 -17.705144528801178 -21.282553458739564 28.804200273643293 -1.201232743965325 1.1816370994722911 -4.6413864468777106 0.16271318451015376 0.20039674491167364 5.0802083002510958 0.34883847783049576 0.79236597310320867 0.13153942668104618 0.024208752100670642];    
    condini = [-38; param(19); param(20); param(21)]
    benchmark=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH11_Cap_Kir_K_L); 
        V=x(1,:);
        benchmark(i,:)=V;
        //plot(t,V)
    end
    //plot(t,benchmark)
endfunction


function [benchmark]=bestFoundRIM_Cap_Kir_K_L_is_Cap_Kir_K_L()
    stimuli=[-15:5:35];
    t0=0;
    t=linspace(0,50,10000);
    
    param= [0.28194480066472999 0.3629941334440131 0.10000000000000092 0.23565162258260272 66.916420143873182 -99.120867483304551 -75.266815972228599 -13.018186452380045 -90 -17.649141234189738 -18.730621009811976 28.555549999949505 -1.2083985275133826 1.9499502620021651 -3.8865297203138178 0.18700000186388682 0.040925646657574283 4.9076914307656727 0.26165548212512507 0.054284485063274796 0.16960424432337129 0.027188874598529693];    
    condini = [-38; param(19); param(20); param(21)]
    benchmark=zeros(11,length(t));
    for i=1:11
        I=stimuli(i);
        x=ode(condini,t0,t,HH11_Cap_Kir_K_L); 
        V=x(1,:);
        benchmark(i,:)=V;
        //plot(t,V)
    end
    //plot(t,benchmark)
endfunction

function plotRIM_Cap_Kir_K_L(plotSolution)
    if ~exists("plotSolution","local") then
         plotSolution = %F;
    end
  
    stimuli=[-15:5:35];
    t=linspace(0,50,10000)
    
    benchmark = benchmarkRIM_Cap_Kir_K_L();
    model     = bestFoundRIM_Cap_Kir_K_L_is_Cap_Kir_K_L();
    
    clf();
    p=gca();
    
    p.foreground=1;
    p.thickness=4;
    p.mark_style=9;
    plot(t,benchmark,'c-');
    f=gcf();
    f.figure_size = [600 730];
    f.children.margins(2)=0.02;

    if plotSolution then
        p.foreground=0;
        //p.thickness=1;
        p.mark_style=9;
        plot(t,model,'k-');
        f=gcf();
        f.children.children(1).children.line_style=3;
        f.children.children(1).children.thickness=2;
        title('$Benchmark: I_{Ca,p}+I_{Kir}+I_K+I_L \quad vs. \quad\\ Model:  I_{Ca,p}+I_{Kir}+I_K+I_L$','fontsize',5)
        //legend(['$Benchmark$';'$Model$'],1)
    else
        title('$Benchmark: \quad I_{Ca,p}+I_{Kir}+I_K+I_L$','fontsize',5)
    end   
        xlabel('$t(10^{-1}s)$','fontsize',5);
        ylabel('$V(t)$','fontsize',5);   
endfunction


