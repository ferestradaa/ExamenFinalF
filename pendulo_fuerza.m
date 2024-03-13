    %Limpieza de pantalla
    
    
    tic 
    %Declaración de variables simbólicas
    syms th1(t)  t  %Angulos de cada articulación
    syms th1p(t) 
    syms th1pp(t) 
    syms m1 m2 Ixx1 Iyy1 Izz1    %Masas y matrices de Inercia
    syms l1 lc1    %l=longitud de eslabones y lc=distancia al centro de masa de cada eslabón
    syms pi g a cero
    
     %Creamos el vector de coordenadas articulares
      Q= [th1];
     %disp('Coordenadas generalizadas');
     %pretty (Q);
     
     %Creamos el vector de velocidades articulares
     Qp= [th1p];
     %disp('Velocidades generalizadas');
     %pretty (Qp);
     %Creamos el vector de aceleraciones articulares
     Qpp= [th1pp]; %se utiliza este formato para simplificar la impresion de resultados
     %disp('Aceleraciones generalizadas');
     %pretty (Qpp);
    
    %Configuración del robot, 0 para junta rotacional, 1 para junta prismática
    RP=[0];
    
    %Número de grado de libertad del robot
    GDL= size(RP,2);
    GDL_str= num2str(GDL);
    
    %Articulación 1 
    %Posición de la articulación 1 respecto a 0
    P(:,:,1)= [l1*cos(th1); l1*sin(th1);0];
    %Matriz de rotación de la junta 1 respecto a 0.... 
    R(:,:,1)= [cos(th1) -sin(th1)  0;
               sin(th1)  cos(th1)  0;
               0         0         1];
    
    
    
    %Creamos un vector de ceros
    Vector_Zeros= zeros(1, 3);
    
    %Inicializamos las matrices de transformación Homogénea locales
    A(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
    %Inicializamos las matrices de transformación Homogénea globales
    T(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
    %Inicializamos las posiciones vistas desde el marco de referencia inercial
    PO(:,:,GDL)= P(:,:,GDL); 
    %Inicializamos las matrices de rotación vistas desde el marco de referencia inercial
    RO(:,:,GDL)= R(:,:,GDL); 
    
    
    for i = 1:GDL
        i_str= num2str(i);
       %disp(strcat('Matriz de Transformación local A', i_str));
        A(:,:,i)=simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);
       %pretty (A(:,:,i));
    
       %Globales
        try
           T(:,:,i)= T(:,:,i-1)*A(:,:,i);
        catch
           T(:,:,i)= A(:,:,i);
        end
    %     disp(strcat('Matriz de Transformación global T', i_str));
        T(:,:,i)= simplify(T(:,:,i));
    %     pretty(T(:,:,i))
    
        RO(:,:,i)= T(1:3,1:3,i);
        PO(:,:,i)= T(1:3,4,i);
        %pretty(RO(:,:,i));
        %pretty(PO(:,:,i));
    end
    
    
    %Calculamos el jacobiano lineal de forma diferencial
    %disp('Jacobiano lineal obtenido de forma diferencial');
    %Derivadas parciales de x respecto a th1 y th2
    Jv11= functionalDerivative(PO(1,1,GDL), th1);
    %Derivadas parciales de y respecto a th1 y th2
    Jv21= functionalDerivative(PO(2,1,GDL), th1);
    %Derivadas parciales de z respecto a th1 y th2
    Jv31= functionalDerivative(PO(3,1,GDL), th1);
    
    %Creamos la matríz del Jacobiano lineal
    %jv_d=simplify([Jv11;
     %             Jv21 ;
     %             Jv31 ]);
    %pretty(jv_d);
    
    %Calculamos el jacobiano lineal de forma analítica
    Jv_a(:,GDL)=PO(:,:,GDL);
    Jw_a(:,GDL)=PO(:,:,GDL);
    
    for k= 1:GDL
        if RP(k)==0 
           %Para las juntas de revolución
            try
                Jv_a(:,k)= cross(RO(:,3,k-1), PO(:,:,GDL)-PO(:,:,k-1));
                Jw_a(:,k)= RO(:,3,k-1);
            catch
                Jv_a(:,k)= cross([0,0,1], PO(:,:,GDL));%Matriz de rotación de 0 con respecto a 0 es la Matriz Identidad, la posición previa tambien será 0
                Jw_a(:,k)=[0,0,1];%Si no hay matriz de rotación previa se obtiene la Matriz identidad
             end
         else
    %         %Para las juntas prismáticas
            try
                Jv_a(:,k)= RO(:,3,k-1);
            catch
                Jv_a(:,k)=[0,0,1];
            end
                Jw_a(:,k)=[0,0,0];
         end
     end    
    
    %Obtenemos SubMatrices de Jacobianos
    Jv_a= simplify (Jv_a);
    Jw_a= simplify (Jw_a);
    %disp('Jacobiano lineal obtenido de forma analítica');
    %pretty (Jv_a);
    %disp('Jacobiano ángular obtenido de forma analítica');
    %pretty (Jw_a);
    
    %Matriz de Jacobiano Completa
    %disp('Matriz de Jacobiano');
    Jac= [Jv_a;
          Jw_a];
    Jacobiano= simplify(Jac);
    %pretty(Jacobiano);
    
    %Obtenemos vectores de Velocidades Lineales y Angulares
    % disp('Velocidad lineal obtenida mediante el Jacobiano lineal');
    V=simplify (Jv_a*Qp);
    % pretty(V);
    % disp('Velocidad angular obtenida mediante el Jacobiano angular');
    W=simplify (Jw_a*Qp);
    %     pretty(W);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Energía Cinética
    
    %Distancia del origen del eslabón a su centro de masa
    %Vectores de posición respecto al centro de masa
     P01=subs(P(:,:,1)/2, l1, lc1); %La función subs sustituye l1 por lc1 en 
                                    %la expresión P(:,:,1)/2
    
    %Creamos matrices de inercia para cada eslabón
    I1=[Ixx1 0 0; 
        0 Iyy1 0; 
        0 0 Izz1];
    
    
    
    %Función de energía cinética
    
    %Extraemos las velocidades lineales en cada eje
    V=V(t);
    Vx= V(1,1);
    Vy= V(2,1);
    Vz= V(3,1);
    
    %Extraemos la velocidad angular en cada ángulo de Euler
    W=W(t);
    W_pitch= W(1,1);
    W_roll= W(2,1);
    W_yaw= W(3,1);
    
    %Calculamos las velocidades para cada eslabón 
    %Ya lo calculamos previamente al multiplicar la matriz jacobiana por Qp
    %Calculamos la energía cinética para cada uno de los eslabones
    
    %Eslabón 1
    %Calculamos el jacobiano lineal y angular de forma analítica
    Jv_a1(:,GDL)=PO(:,:,GDL);
    Jw_a1(:,GDL)=PO(:,:,GDL);
    
    for k= 1:GDL
        if RP(k)==0 
           %Para las juntas de revolución
            try
                Jv_a1(:,k)= cross(RO(:,3,k-1), PO(:,:,GDL)-PO(:,:,k-1));
                Jw_a1(:,k)= RO(:,3,k-1);
            catch
                Jv_a1(:,k)= cross([0,0,1], PO(:,:,GDL));%Matriz de rotación de 0 con respecto a 0 es la Matriz Identidad, la posición previa tambien será 0
                Jw_a1(:,k)=[0,0,1];%Si no hay matriz de rotación previa se obtiene la Matriz identidad
             end
         else
    %         %Para las juntas prismáticas
            try
                Jv_a1(:,k)= RO(:,3,k-1);
            catch
                Jv_a1(:,k)=[0,0,1];
            end
                Jw_a1(:,k)=[0,0,0];
         end
     end    
    
    %Obtenemos SubMatrices de Jacobianos
    Jv_a1= simplify (Jv_a1);
    Jw_a1= simplify (Jw_a1);
    %disp('Jacobiano lineal obtenido de forma analítica');
    %pretty (Jv_a);
    %disp('Jacobiano ángular obtenido de forma analítica');
    %pretty (Jw_a);
    
    %Matriz de Jacobiano Completa
    %disp('Matriz de Jacobiano');
    Jac1= [Jv_a1;
          Jw_a1];
    Jacobiano1= simplify(Jac1);
    % pretty(Jacobiano);
    
    %Obtenemos vectores de Velocidades Lineales y Angulares
     %disp('Velocidad lineal obtenida mediante el Jacobiano lineal del Eslabón 1');
    Qp=Qp(t);
    V1=simplify (Jv_a1*Qp(1));
     %pretty(V1);
     % disp('Velocidad angular obtenida mediante el Jacobiano angular del Eslabón 1');
    W1=simplify (Jw_a1*Qp(1));
     % pretty(W1);
    
    %Calculamos la energía cinética para cada uno de los eslabones%%%%%%%%%%
    
    %Eslabon1
    V1_Total= V1+cross(W1,P01);
    K1= (1/2*m1*(V1_Total))'*((V1_Total)) + (1/2*W1)'*(I1*W1);
    disp('Energía Cinética en el Eslabón 1');
    K1= simplify (K1);
    pretty (K1)
    
    
    
    K_Total= simplify (K1)
    pretty (K_Total)
    %K_Total= simplify (K1+K2);
    
    %Calculamos la energía potencial para cada uno de los eslabones
    
    %Obtenemos las alturas respecto a la gravedad
     h1= P01(2); %Tomo la altura paralela al eje z
     
     U1=m1*g*h1;
    
    
     %Calculamos la energía potencial total
     U_Total= U1 
    
     %Obtenemos el Lagrangiano
     Lagrangiano= simplify (K_Total-U_Total)
     pretty (Lagrangiano);
    
    %Modelo de Energía
     H= simplify (K_Total+U_Total);
      pretty (H)
    
    
    
      %%%%%%%%%%%%%%%%%%%%%%Ecuaciones de Movimiento%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     %Lagrangiano derivado con respecto a la primera coordenada generalizada de
     %velocidad
    
     %Definimos un vector columna de derivadas con respecto al tiempo
     %En este vector agrego las velocidades y aceleraciones
      %Derivadas respecto al tiempo
     Qd=[th1p(t); th1pp(t)];
    
     %Obtenemos las derivadas de la velocidad en la primera coordenada
     %generalizada
     dQ1=[diff(diff(Lagrangiano,th1p), th1),... %Derivamos con respecto a la primera velocidad generalizada th1p para las 3 posiciones articulaciones
         diff(diff(Lagrangiano,th1p), th1p)];%Derivamos con respecto a la primera velocidad generalizada th1p para las 3 velocidades articulaciones
    
     %Definimos el torque 1
     t1= dQ1*Qd- diff(Lagrangiano, th1);
    
     %Generación del Modelo Dinámico en forma matricial
    
    %Matriz de Inercia
    
    %Extraemos coeficientes de aceleraciones
    
    M=[diff(t1, th1pp)];
    rank (M);
    
    M=M(t);
    
    
    
     %Fuerzas Centrípetas y de Coriolis
     
     %Definimos Mp
    
     Mp=[diff(M(1,1),th1)]*Qp;%Se deriva parcialmente en el tiempo respecto a todas las variables 
    
    %Definimos la energía cinética en su forma matricial
    k=1/2*transpose(Qp)*M*Qp;
    
    %Definimos dk
    
    dk=[diff(k, th1)];
    
    %Fuerzas centrípetas y de Coriolis
     C= Mp*Qp-dk;
    
    
     %Par Gravitacional
     %se sustituyen las velocidades y acele raciones por 0
     r=cero;
     a1=subs(t1, th1p, r);
    
     %Torque gravitacional en el motor 1
     G1=a1;
    
    % Vector de par gravitacional
    
    G=[G1];
    
    toc

