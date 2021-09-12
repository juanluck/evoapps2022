/////////////////////////////////////////////////////////////
///////////////     Récupération données      ///////////////
/////////////////////////////////////////////////////////////

t=linspace(0,50,10000);
t0=0;
stimulus=[-15:5:35];

/////////////////////////////////////////////////////////////
///////////////    Fonction coût algorithme    //////////////
/////////////////////////////////////////////////////////////


function y=costfunction(modelVal,benchmarkVal)
    c=0;
    
    traces=size(modelVal)(1);
    points=size(modelVal)(2);    
    for i=1:traces
        for k=1:points
            c=c+(modelVal(i,k)-benchmarkVal(i,k))*(modelVal(i,k)-benchmarkVal(i,k));
        end
    end
    y=sqrt(c)/points;
endfunction


////////////////////////////////////////////////////////
/////////    Estimation de la capacitance C    /////////
////////////////////////////////////////////////////////

function [bM, valBest]=optimization(model,constraints,benchmark)
    NP=140;
    itermax=1000;
    F=0.5;
    CR=0.9;
    
    //disp("Step 1")
    [benchmarkVal]=benchmark();
    //disp("OK")    
    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////
    [Xmin,Xmax]=constraints();
    D=length(Xmin);

    pop=zeros(D,NP);
    
    /////////////////////////////////////////
    //// Initialisation de ma population ////
    /////////////////////////////////////////

    for j=1:NP
        for i=1:D
            pop(i,j)=Xmin(i)+(Xmax(i)-Xmin(i))*rand();
        end
    end
    
    
    //////////////////////////////////////////////////////////////
    //// Évaluation du meilleur individu après initialisation ////
    //////////////////////////////////////////////////////////////
    
    val=zeros(NP,1); // tableau avec le coût de chacun des individus
    //disp("Step 2")
    for j=1:NP
        [modelVal,pop(:,j)]=model(pop(:,j));
        val(j)=costfunction(modelVal,benchmarkVal);
    end
    //disp("OK")

    

    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
    
    //disp("Init pop")
    //disp(val(bestIndex))
    
    
    
    ////////////////////////
    //// Étape suivante ////
    ////////////////////////
     
    iter=1; // nombre d'itération
    U=zeros(D,NP); // Vecteur intermédiaire perturbé (mutation + crossover)
    tempval=0;
    
    //timestamp=getdate("s");
    while iter<itermax
        for j=1:NP
            // ======= Construction de la matrice U = variation différentielle + crossover =======

            // ========= Tirage aléatoire de 3 entiers distincts r1, r2 et r3 et différents de j ========
            r1=j; r2=j; r3=j;//////////////////////////////////////
            while (r1==r2 | r1==r3 | r2==r3 | r1==j | r2==j | r3==j)
                r1=floor(1+NP*rand());
                r2=floor(1+NP*rand());
                r3=floor(1+NP*rand());
            end
            // ======== Variation différentielle =======
            V=pop(:,r1) + F*(pop(:,r2)-pop(:,r3));
            
            for i=1:D
                if V(i)<=Xmin(i) then V(i)=Xmin(i);
                elseif V(i)>Xmax(i) then V(i)=Xmax(i);
                end
            end
            // ======== Crossover ========
            for i=1:D
                if rand()<CR then
                    U(i,j)=V(i);
                else
                    U(i,j)=pop(i,j);
                end
            end
        end // fin for j=1:NP
    
    // ======== Sélection ========
        for j=1:NP
            [modelVal,U(:,j)]=model(U(:,j));
            tempval=costfunction(modelVal,benchmarkVal);

            if tempval<=val(j) then
                pop(:,j) = U(:,j);
                val(j) = tempval;
            end
        end
        //disp(iter)
        iter = iter + 1;
        bestIndex=1;
        for b=2:NP
            if val(b)<val(bestIndex) then bestIndex=b; end
        end
        costVec(iter)=val(bestIndex);
        //disp(costVec(iter));
    end  //fin de la boucle while
    //timestamp=(getdate("s")-timestamp);
    
    // Détermination de l'indice du meilleur individu
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
//    disp(bestIndex);
    
    // Sauvegarde du meilleur individu
    bM = [];
    bM = pop(:,bestIndex);
    
    //disp(val);
    //disp(bM);
    //disp(val(bestIndex));
    valBest=val(bestIndex);
    //disp("Seconds-->",timestamp)
//    iterVec=1:1:itermax;
//    plot(iterVec,costVec,2)
endfunction






