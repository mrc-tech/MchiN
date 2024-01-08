% Questa classe è basata sul codice SeReS_strisce

classdef MchiNstrisceSteel < handle
    properties
        hmax %altezza massima del layer
        layer %matrice dove ogni riga è un layer e le colonne sono: y, h, b_bot, b_top
%         materials %struttura che contiene le funzioni sigma(eps) dei materiali (ogni riga è un materiale differente)
        legameAcc %legame sigma(eps) dell'acciaio
        %proprietà aggiuntive:
        layersProp %ogni riga è un layer, le colonne sono: Area, ygl
        %proprietà dei materiali:
        fy %resistenza acciaio
        epssy, epssu %deformazione di snervamento e utlima dell'acciaio
        %costanti della sezione:
        Area    %Area della Sezione
        yg      %coordinata del baricentro
        Ix      %Momento di inerzia secondo l'asse orizzontale
        H       %altezza totale della sezione
    end
    methods
        function obj = MchiNstrisceSteel(hmax,Acc)
            %costruttore della classe
            obj.hmax = hmax;
            obj.legameAcc = Acc;
        end
        
        function addLayer(obj,h,b,t)
            %aggiunge tanti layer alti massimo h_max fino a creare la patch
            %desiderata
            y = 0;
            for i=1:size(obj.layer,1); y = y + obj.layer(i,2); end
            sum = 0;
            while (h - (sum+obj.hmax) > 0)
                bj = b + (t-b)/(h)*(sum);
                tj = b + (t-b)/(h)*(sum+obj.hmax);
                if size(obj.layer,1) > 0; obj.layer = [obj.layer; (obj.layer(size(obj.layer,1),1)+obj.layer(size(obj.layer,1),2)) obj.hmax bj tj];
                else; obj.layer = [obj.layer; 0 obj.hmax b tj]; end %aggiunge il primo layer
                sum = sum + obj.hmax;
            end
            bj = b + (t-b)/(h)*(sum);
            if size(obj.layer,1) > 0; obj.layer = [obj.layer; (obj.layer(size(obj.layer,1),1)+obj.layer(size(obj.layer,1),2)) h-sum bj t];
            else; obj.layer = [obj.layer; 0 h b t]; end %aggiunge il primo layer
        end
        
        function defineMaterial(obj,fy,epssy,epssu)
            %aggiunge le proprietà dei materiali
            if nargin < 4
                obj.fy = fy;
                obj.epssy = 0.002;
                obj.epssu = 0.15;
            else
                obj.fy = fy;
                obj.epssy = epssy;
                obj.epssu = epssu;
            end
        end
        
        function initgeo(obj)
            %calcola la geometria della sezione
            %si prende in considerazione il fatto che sia composta tutta dello stesso materiale (in realtà dovrei SEMPRE omogeneizzare i layer [per ogni layer un coeff. di omogenizzazione n])
            obj.Area = 0;
            obj.Ix = 0;
            Sx = 0;
            for i=1:size(obj.layer,1)
                %1  2  3  4 
                %y, h, b, t
                Areal = 0.5 * obj.layer(i,2) * (obj.layer(i,3) + obj.layer(i,4));
                ygl   = obj.layer(i,1) + obj.layer(i,2)/3 * (1+obj.layer(i,4)/(obj.layer(i,3)+obj.layer(i,4)));
                obj.Area = obj.Area + Areal;
                obj.layersProp = [obj.layersProp; Areal ygl];
                Sx = Sx + obj.layer(i,3) * obj.layer(i,2) * (obj.layer(i,2)/2 + obj.layer(i,1)) + 0.5*obj.layer(i,2)*(obj.layer(i,4)-obj.layer(i,3))*(2/3.*obj.layer(i,2) + obj.layer(i,1));
                obj.Ix = obj.Ix + 1/36*(obj.layer(i,2))^3*(2*obj.layer(i,3)+obj.layer(i,4)) + obj.layer(i,2)*(obj.layer(i,3)*(1/2*obj.layer(i,2)+obj.layer(i,1))^2 + 1/2*(obj.layer(i,4)-obj.layer(i,3))*(2/3*obj.layer(i,2)+obj.layer(i,1))^2);
            end
            obj.yg = Sx / obj.Area;
            obj.Ix = obj.Ix - obj.Area * (obj.yg)^2;
            obj.H = sum(obj.layer(:,2));
        end
        
        function [N,M] = setStrain(obj,eps_b,eps_t)
            %calcola N e M impostata ladeformazione al lembo inferiore
            %(spe_b) e al lembo superiore (eps_t)
            
            %forze nei layers:
            F = obj.legameAcc((eps_t-eps_b)*obj.layersProp(:,2)/obj.H+eps_b, obj.epssy,obj.epssu,obj.fy) .* obj.layersProp(:,1);
            
            %contributo dell'acciaio
            N = sum(F);
            M = sum(F .* obj.layersProp(:,2));
            
            M = M - N * obj.yg; %aggiunge il momento dovuto alla forza normale
        end
        
        function [M,phi,yn] = yieldingPoint(obj, N)
            %calcola il punto di snervamento dato uno sforzo normale N
            
            if N < obj.fc*obj.Area; error('la sezione non può resistere a tale compressione'); end
            if N > obj.fy*sum(obj.rebar(:,2)); error('la sezione non può resistere a tale trazione'); end
            
            % fissato lo eps_inf = epssy cambia eps_top fino a trovare N corretto
            eps_inf = obj.epssy;
            
            %trova un range per fare la bisezione
            count = 0;
            eps_supP = obj.epssy; %l'ho scelto io (l'N iniziale può essere solamente minore)
            eps_supN = obj.epssy; %valore Negativo
            while obj.setStrain(eps_inf,eps_supN) > N
                eps_supN = eps_supN - 0.001;
                count = count + 1;
                if count > 1000; error('MchiNstrisce non è arrivato a convergenza'); end
            end
            
            %fa la bisezione
            err = [100 100 100];
            count = 0;
            while ((abs(err(2)) > abs(0.005 * N)) )%|| (abs(eps_supP-eps_supN) > abs(0.0001*eps_supN))) %+0.000001 serve per non creare errori per N=0
                err = [(obj.setStrain(eps_inf,eps_supP) - N) (obj.setStrain(eps_inf,(eps_supP+eps_supN)/2) - N) (obj.setStrain(eps_inf,eps_supN) - N)];
                if err(1)*err(2) < 0
                    %la risultante è tra err(1) e err(2) (eps_supP è corretto)
                    eps_supN = (eps_supP+eps_supN)/2;
                    
                else
                    %la risultante è tra err(2) e err(3) (eps_supN è corretto)
                    eps_supP = (eps_supP+eps_supN)/2;
                end
                count = count + 1;
                if count > 1000; error('MchiNstrisce non è arrivato a convergenza'); end
            end
            
            eps_sup = (eps_supP+eps_supN)/2;
            
            [~,M] = obj.setStrain(eps_inf,eps_sup);
            phi = (eps_inf-eps_sup)/obj.H; %compressione: negativa (la fibra superiore è compressa)
            yn = -eps_sup/(eps_inf-eps_sup) * obj.H;
            M = -M; %convenzione italiana
        end
        
%         function [M,phi,yn] = ultimatePoint(obj, N)
%             %calcola il momento e la curvatura ultimi
%             
%             if N < obj.fc*obj.Area; error('la sezione non può resistere a tale compressione'); end
%             if N > obj.fy*sum(obj.rebar(:,2)); error('la sezione non può resistere a tale trazione'); end
%             
%             % fissato lo eps_top = epscu cambia eps_top fino a trovare N corretto
%             eps_sup = obj.epscu;
%             
%             %trova un range per fare la bisezione
%             count = 0;
%             eps_infP = obj.epscu; %l'ho scelto io sta sicuramente qua dentro...
%             eps_infN = obj.epscu; %in questo intevallo (non può essere più compresso di così)
%             rottura_acciaio = 0;
%             while obj.setStrain(eps_sup,eps_infP) < N
%                 Nint = obj.setStrain(eps_sup,eps_infP);
%                 eps_infP = eps_infP + 0.001;
%                 %l'acciaio è arrivato a rottura, quindi si ha una rottura lato acciaio
%                 if eps_infP > obj.epssu
%                     eps_infP = obj.epssu;
%                     eps_sup = eps_sup + 0.0003;
%                     rottura_acciaio = 1;
%                 end
%                 count = count + 1;
%                 if count > 1000; error('MchiNstrisce non è arrivato a convergenza'); end
%             end
% %             if rottura_acciaio; warning('rottura lato acciaio'); end
%             if rottura_acciaio; logFile('MchiN::ultimatePoint: rottura lato acciaio'); end
%             
%             %fa la bisezione
%             err = [100 100 100];
%             count = 0;
%             while ((abs(err(2)) > abs(0.000001 * N)) || (abs(eps_infP-eps_infN) > abs(0.0001*eps_infN))) %+0.000001 serve per non creare errori per N=0)
%                 err = [(obj.setStrain(eps_infP,eps_sup) - N) (obj.setStrain((eps_infP+eps_infN)/2,eps_sup) - N) (obj.setStrain(eps_infN,eps_sup) - N)];
%                 if err(1)*err(2) < 0
%                     %la risultante è tra err(1) e err(2) (eps_supP è corretto)
%                     eps_infN = (eps_infP+eps_infN)/2;
%                     
%                 else
%                     %la risultante è tra err(2) e err(3) (eps_supN è corretto)
%                     eps_infP = (eps_infP+eps_infN)/2;
%                 end
%                 count = count + 1;
%                 if count > 10000; error('MchiNstrisce non è arrivato a convergenza'); end
%                 %controlli
%                 [Np,~] = obj.setStrain(eps_infP,eps_sup);
% %                 [Nn,~] = obj.setStrain(eps_infN,eps_sup);
% %                 [Nm,~] = obj.setStrain((eps_infN+eps_infP)/2,eps_sup);
%                 if Np < N
% %                     disp('qualcosa non va proprio');
%                     eps_sup = eps_sup + 0.0001;
%                     eps_infP = obj.epssu;
%                 end
%                 if count > 50
%                     disp(['sta succendendo qualcosa di brutto' num2str(N) ' ' num2str(Np) ' ' num2str(err(2))]);
%                     eps_sup = eps_sup - 0.0001; %lo ha aumentato troppo...
%                     if eps_sup < obj.epscu
%                         %resetta tutto
%                         eps_sup = obj.epscu;
%                         eps_infN = obj.epscu;
%                         eps_infP = obj.epssu;
%                     end
%                 end
%                 if count > 80
%                     disp('STA PER FALLIRE:: ABORTING...');
%                     err(2) = 0; %abortisce il calcolo
%                 end
%             end
%             
%             eps_inf = (eps_infP+eps_infN)/2;
%             
%             [~,M] = obj.setStrain(eps_inf,eps_sup);
%             phi = (eps_inf-eps_sup)/obj.H; %compressione: negativa (la fibra superiore è compressa)
%             yn = -eps_sup/(eps_inf-eps_sup) * obj.H;
%             M = -M; %convenzione italiana
%         end
        
        function [M,phi,yn] = ultimatePoint(obj, N)
            %calcola il momento e la curvatura ultimi
            
            if abs(N) > obj.fy*obj.Area; error('la sezione non può resistere a tale sforzo'); end
            
            % fissato lo eps_top = epscu cambia eps_top fino a trovare N corretto (se la rottura è lato CLS)
            eps_sup = obj.epscu;
            
            %trova un range per fare la bisezione:
            count = 0;
            eps_infP = obj.epscu; %l'ho scelto io sta sicuramente qua dentro...
            eps_infN = obj.epscu; %in questo intevallo (non può essere più compresso di così)
            eps_supP = obj.epscu; %nel caso si ha la rottura lato acciaio la bisezione va fatta sulla eps_sup
            eps_supN = obj.epscu;
            rottura_acciaio = 0;
            while obj.setStrain(eps_infP,eps_supP) < N
                eps_infP = eps_infP + 0.0001;
                %l'acciaio è arrivato a rottura, quindi si ha una rottura lato acciaio
                if eps_infP > obj.epssu
                    eps_infP = obj.epssu; %fissa eps_inf
                    eps_inf = eps_infP; %usa un solo eps_inf (fissato) per non confondere le idee...
                    eps_supP = eps_supP + 0.001;
                    rottura_acciaio = 1;
                end
                count = count + 1;
                if count > 1000; error('MchiNstrisce non è arrivato a convergenza'); end
            end
%             if rottura_acciaio; warning('rottura lato acciaio'); end
            if rottura_acciaio; logFile('MchiN::ultimatePoint: rottura lato acciaio'); end
            
            %fa la bisezione
            err = [100 100 100];
            count = 0;
            if rottura_acciaio == 0
                %rottura del CLS:
                while ((abs(err(2)) > abs(0.000001 * N)) || (abs(eps_infP-eps_infN) > abs(0.0001*eps_infN))) %+0.000001 serve per non creare errori per N=0)
                    err = [(obj.setStrain(eps_infP,eps_sup) - N) (obj.setStrain((eps_infP+eps_infN)/2,eps_sup) - N) (obj.setStrain(eps_infN,eps_sup) - N)];
                    if err(1)*err(2) < 0
                        %la risultante è tra err(1) e err(2) (eps_infP è corretto)
                        eps_infN = (eps_infP+eps_infN)/2;
                    else
                        %la risultante è tra err(2) e err(3) (eps_infN è corretto)
                        eps_infP = (eps_infP+eps_infN)/2;
                    end
                    count = count + 1;
                    if count > 10000; error('MchiNstrisce non è arrivato a convergenza'); end
                end
                eps_inf = (eps_infP+eps_infN)/2;
            else
                %rottura dell'acciaio:
                while ((abs(err(2)) > abs(0.000001 * N)) || (abs(eps_supP-eps_supN) > abs(0.0001*eps_supN))) %+0.000001 serve per non creare errori per N=0)
                    err = [(obj.setStrain(eps_inf,eps_supP) - N) (obj.setStrain(eps_inf,(eps_supP+eps_supN)/2) - N) (obj.setStrain(eps_inf,eps_supN) - N)];
                    if err(1)*err(2) < 0
                        %la risultante è tra err(1) e err(2) (eps_supP è corretto)
                        eps_supN = (eps_supP+eps_supN)/2;
                    else
                        %la risultante è tra err(2) e err(3) (eps_supN è corretto)
                        eps_supP = (eps_supP+eps_supN)/2;
                    end
                    count = count + 1;
                    if count > 9000
                        disp('sta per fallire');
                    end
                    if count > 10000; error('MchiNstrisce non è arrivato a convergenza'); end
                end
                eps_sup = (eps_supP+eps_supN)/2;
            end
            
            [~,M] = obj.setStrain(eps_inf,eps_sup);
            phi = (eps_inf-eps_sup)/obj.H; %compressione: negativa (la fibra superiore è compressa)
            yn = -eps_sup/(eps_inf-eps_sup) * obj.H;
            M = -M; %convenzione italiana
        end
        
        function eps = deformazione(obj,eps0,chi,y)
            %funzione di epsilon: def(eps0,chi,y) = eps0 + chi*(y-yg)
            eps = eps0 + chi*(y-obj.yg);
        end
        
        function [f] = F(obj,eps0,chi,N0,M0)
            %calcola le funzioni
            %FUNZIONE AUSILIARIA
            [N,M]  = obj.Fint(eps0,chi);
            f(1,1) = N - N0;
            f(2,1) = M - M0;
        end
        
        function [N,M] = Fint(obj,eps0,chi)
            %restituisce le forze interne della sezione dato un profilo di
            %deformazione
            
            N = sum(obj.legameAcc(obj.deformazione(eps0,chi,(obj.layersProp(:,2))), obj.epssy,obj.epssu,obj.fy) .* obj.layersProp(:,1));
            M = sum(obj.legameAcc(obj.deformazione(eps0,chi,(obj.layersProp(:,2))), obj.epssy,obj.epssu,obj.fy) .* (obj.layersProp(:,2)-obj.yg) .* obj.layersProp(:,1));
        end
        
        function [eps0,chi] = findEquilibrium(obj,N,M)
            %trova il profilo di deformazione che equilibria le azione
            %esterne N e M. Ha l'ipotesi di sezion piane che rimangono
            %piane anche dopo la deformazione quindi eps0 è la deformazione
            %in corrispondenza della fibra baricentrica e chi è la
            %curvatura.
            %usa il metodo di Newton-Raphson
            
            Phi= [0; 0]; %vettore eps,chi
            err = 1; %errore che deve andare a zero
            tol = 1.e-10; %tolleranza
            delta = 1.e-8;
            
            while err > tol
                %calcola il jacobiano
                JF(:,1) = (obj.F(Phi(1)+delta,Phi(2),N,M)-obj.F(Phi(1),Phi(2),N,M))/delta;
                JF(:,2) = (obj.F(Phi(1),Phi(2)+delta,N,M)-obj.F(Phi(1),Phi(2),N,M))/delta;
                %calcola lo step successivo
                Phi0=Phi;
                Phi = Phi0 - inv(JF) * (obj.F(Phi0(1),Phi0(2),N,M));
                %calcola l'errore
                err = norm(Phi - Phi0)/norm(Phi0);
            end
            
%             [~,M] = obj.setStrain(Phi(1)-Phi(2)*obj.yg, Phi(1)+Phi(2)*(obj.H-obj.yg));
            eps0 = Phi(1);
            chi  = Phi(2);        
            
        end
        
%         function [M,chi] = curvaMchi(obj,N)
%             %calcola la curva del momento in funzione della curvatura
%             %usa il metodo di Newton
%             
%             M = linspace(0,0.05,100);
%             
%             for i=1:size(M,2)
%                 [eps0,chi(i)] = obj.findEquilibrium(N,M(i));
%             end
%             
%         end
        
        function [M,chi] = curvaMchi2(obj,N,chi_lim)
            %calcola la curva M-chi-N
            
            if nargin < 3; chi_max = (2*obj.epssu)/obj.H; else; chi_max = chi_lim; end
            
            N_punti = 100; %numero punti della curva
            delta = 1.e-8; %variazione per il calcolo numerico della derivata parziale
            tol = 1.e-3; %tolleranza per la convergenza
            
            chi = linspace(0,chi_max,N_punti);
            M   = zeros(1,N_punti);
            
            for i=2:(N_punti) %parte da 2 perchè così M(0) = 0 (altrimenti non lo è per problemi numerici...)
                err   = 1;
                eps0  = 0;
                count = 0;
                bisez = 0; %controlla se non bisogna ricorrere alla bisezione
                Nint = 0;
                deltaNint = 100;
                while (err > tol) || (deltaNint > 0.001*abs(N))
                    Nint_old = Nint;
                    [Nint,~] = obj.Fint(eps0,chi(i));
                    deltaNint = abs(Nint-Nint_old);
                    derivata = (obj.Fint(eps0+delta,chi(i))-obj.Fint(eps0,chi(i)))/delta;
                    eps0_old = eps0;
                    eps0 = eps0 + (1/derivata) * (N - Nint);
                    if eps0 ~= 0; err = abs((eps0 - eps0_old)/eps0); else; err = abs(eps0 - eps0_old); end
                    count = count + 1;
                    if count > 100 || derivata == 0
%                         error('MchiN non è arrivato a convergenza');
                        logFile('MchiN::curvaMchi2: Newton-Raphson non va bene, calcolo con la bisezione...');
                        bisez = 1;
                        break;
                    end
                end
                if bisez
                    %cerca il range di laboro per la bisezione:
                    eps0_min = -obj.epssu; %non può essere minore di questo
                    eps0_max = -obj.epssu;
                    Nint = obj.Fint(eps0_max,chi(i)); %ATTENZIONE: se la curvatura è elevata questo potrebbe essere maggiore di N (in tutti gli altri casi DOVREBBE essere MINORE)
                    if Nint < N
                        while Nint < N
                            eps0_max = eps0_max + 0.0001;
                            Nint = obj.Fint(eps0_max,chi(i));
                        end
                    else
                        while Nint > N
                            eps0_max = eps0_max + 0.000001;
                            Nint = obj.Fint(eps0_max,chi(i));
                            if eps0_max > obj.epssu
                                %C'è QUALCOSA CHE NON VA PROPRIO
                                disp('ERRORE nella ricerca del range per la bisezione (per Nint > N)');
                                logFile('MchiN::curvaMchi2::bisezione: ERRORE nella ricerca del range per la bisezione (per Nint > N)');
                                M(i) = 0;
                                return;
                            end
                        end
                    end
%                     eps0_max = eps0_max + 0.0001; %aumenta eps0_max per motivi numerici (così gli estremi dell'intervallo sono più distanti dal punto)
                    %esegue la bisezione:
                    count = 0;
                    Ntot = [1 1 1];
                    tol_bis = 0.001; %tolleranza per la bisezione
                    while (norm(Ntot) > (tol_bis * abs(N))) || (abs(abs(eps0_max)-abs(eps0_min)) > (tol_bis * abs(eps0_min)))
                        eps0 = (eps0_min + eps0_max) / 2; %valore medio
                        Ntot(1) = obj.Fint(eps0_min,chi(i)) - N;
                        Ntot(2) = obj.Fint(eps0,chi(i)) - N;
                        Ntot(3) = obj.Fint(eps0_max,chi(i)) - N;
                        if (Ntot(1)*Ntot(2)) < 0
                            eps0_max = eps0;
                        else
                            eps0_min = eps0;
                        end
                        count = count + 1;
                        if count > 100
                            disp('ERRORE'); logFile('MchiN::curvaMchi2::bisezione: ERRORE convergenza');
%                             error('MchiN non è arrivato a convergenza');
                            eps0_max = eps0_min; %si INVENTA un risultato (tanto un punto solo toppato non influenza molto la ricerca di chi_y tramite minimi quadrati)
                            Ntot = 0; %USCITA FORZATA
                        end
                    end
                end
                [~,M(i)] = obj.Fint(eps0,chi(i));
            end
            
        end
        
        function [My,chiy,Mu,chiu] = findPoints(obj,N)
            %trova i punti di snervamento e ultimo della sezione dalla
            %curva curvatura-momento.
            
            [~,chi_max] = obj.ultimatePoint(N); %METODO PER STABILIZZARE I RISULTATI NUMERICI
            [M,chi] = obj.curvaMchi2(N,chi_max); %calcola la curva "grezza"
            
            %trova il punto ultimo come il punto in cui la curva comincia
            %ad avere la derivata negativa (i materiali non hanno tratti
            %degradanti a tangente negativa):
%             i = find(diff(M) < 0 ,1); %tronca quando comincia a decrescere CERTE VOLTE FALLISCE PER PROBLEMI NUMERICI ALL'INIZIO
            [~,i] = max(M); %tronca quando raggiunge il massimo
            Mu   = M(i);
            chiu = chi(i);
            %elimina la parte di curva dopo il punto ultimo:
            M   = M(1:i);
            chi = chi(1:i);
            
            %trova il punto di snervamento bilinearizzando la curva:
            P0 = [chiu/3 3/4*Mu]; %parametri iniziali di prova
            model = @(P,x) (P(2)/P(1)*x).*heaviside(P(1)-x) + ((Mu-P(2))/(chiu-P(1))*(x-P(1))+P(2)).*heaviside(x-P(1)); %crea il modello da fittare
            lb = [0 0]; %lower bound (non ci possono essere valori negativi)
            ub = [chiu Mu]; %upper bound (il punto di snervamento non può essere maggiore di quello ultimo)
            options = optimoptions('lsqcurvefit','Display','off');
            Pfit = lsqcurvefit(model,P0,chi,M,lb,ub,options); %fitta il modello minimizzando i quadrati
            chiy = Pfit(1);
            My   = Pfit(2);
            
            %controlla i risultati:
            
            
            
            %comandi per i test:
%             plot(chi,M); hold on;
%             plot(linspace(0,chiu,100),model(Pfit,linspace(0,chiu,100)),'k--');
%             plot(chiy,My,'bo');
%             plot(chiu,Mu,'ro');
            
            
        end
        
        function [My,chiy,Mu,chiu] = findPoints2(obj, N)
            %trova i punti trovando prima chi_u e poi discretizzando la
            %curva in un numero fisso di punti
            
            [Mu,chiu] = obj.ultimatePoint(N);
            [M,chi] = obj.curvaMchi2(N,chiu); %calcola la curva
            
            %trova il punto di snervamento bilinearizzando la curva:
            P0 = [chiu/3 3/4*Mu]; %parametri iniziali di prova
            model = @(P,x) (P(2)/P(1)*x).*heaviside(P(1)-x) + ((Mu-P(2))/(chiu-P(1))*(x-P(1))+P(2)).*heaviside(x-P(1)); %crea il modello da fittare
            lb = [0 0]; %lower bound (non ci possono essere valori negativi)
            ub = [chiu Mu]; %upper bound (il punto di snervamento non può essere maggiore di quello ultimo)
            options = optimoptions('lsqcurvefit','Display','off');
            Pfit = lsqcurvefit(model,P0,chi,M,lb,ub,options); %fitta il modello minimizzando i quadrati
            chiy = Pfit(1);
            My   = Pfit(2);
            
            %DEBUGGING
%             plot(chi,M); hold on;
%             plot(chiy,My,'ob');
%             plot(chiu,Mu,'or');
            
        end
        
        function out = test(obj,chi)
            %esegue dei test sulla sezione
            
            eps0 = -0.003;
%             chi = 0.0001;
            
            for i=1:60
                eps0=eps0+0.0001;
                [N(i),M(i)] = obj.Fint(eps0,chi);
            end
            
            out(:,1) = N;
            out(:,2) = M;
            
        end
        
        function initRect(obj,h,b,As1,As2,c)
            %inizializza una sezione rettangolare
            if nargin < 6; c = 0.03; end
            obj.addLayer(h,b,b); %crea la sezione in CLS
            obj.addRebar(c,As1);
            obj.addRebar(h-c,As2);
            
            obj.initgeo;
        end
        
        function initRect2(obj,h,b,rhotot,rhoc_rhotot,db,c)
            %inizializza una sezione rettangolare (con anche l'armatura di parete)
            if nargin < 7; c = 0.03; end
            Ac = b*h; %area calcestruzzo
            Ab = pi/4 * db^2; %area della barra
            np = max(round((h-2*c)/0.3-1), 0); %calcola il numero di armature di parete (non ci può essere in interferro maggiore di 30 cm tra esse) [maggiore uguale a zero..]
            n2 = round(Ac*rhotot/Ab * rhoc_rhotot); %calcola il numero di armature superiori (compresse)            
            n1 = round(Ac*rhotot/Ab - n2 - 2*np);  %calcola il numero di armature inferiori (tese)
            
            obj.addLayer(h,b,b); %crea la sezione in CLS
            obj.addRebar(c,n1*Ab); %crea l'armatura inferiore
            obj.addRebar(h-c,n2*Ab); %crea l'armatura superiore
            %aggiunge le armature di parete:
            if np > 0
                interferro_laterale = (h-2*c)/(np+1);
                y = c + interferro_laterale;
                for i=1:np
                    obj.addRebar(y,2*Ab);
                    y = y + interferro_laterale;
                end
            end
            
            obj.initgeo; %calcola la geometria della sezione
        end
        
        function initCirc(obj,D,db,rhotot,c)
            %inizializza una sezione circolare
            if nargin < 6; c = 0.03; end
            AreaBarra = pi/4 * db^2;
            numeroBarre = round(rhotot * (D^2)/(db^2));
            
            hs = D/50;
            
            for y=linspace(0,D-hs,D/hs)
                thetab = acos(1 - 2/D*y);
                b = D*sin(thetab);
                thetat = acos(1 - 2/D*(y+hs));
                t = D*sin(thetat);
                obj.addLayer(hs,real(b),real(t));
            end
            for theta=linspace(0,pi,(numeroBarre+1)/2)
                y = D/2 - ((D-2*c)/2)*cos(theta);
                obj.addRebar(y,2*AreaBarra);
            end
            
            obj.initgeo;
        end
        
        function initHollow(obj,h,b,s,db,rhotot,c)
            %inizializza una sezione rettangolare cava
            if nargin < 7; c = 0.03; end
            AreaBarra = pi/4 * db^2;
            nBarre = round(rhotot*(8*s*(h+b-2*s))/(pi*db^2));
            interferro = 4*(h+b-2*s)/nBarre;
            nBarreParete = round(4*(h-2*s)/interferro); %round(8/pi * ((h-2*s)*s)/(db^2));
            nBarreStrati = round((nBarre-nBarreParete)/4);
            
            %crea la sezione in CLS
            obj.addLayer(s,b,b);
            obj.addLayer(h-2*s,2*s,2*s);
            obj.addLayer(s,b,b);
            %aggiunge gli strati orizzontali
            obj.addRebar(c,AreaBarra*nBarreStrati);
            obj.addRebar(s-c,AreaBarra*nBarreStrati);
            obj.addRebar(h-s+c,AreaBarra*nBarreStrati);
            obj.addRebar(h-c,AreaBarra*nBarreStrati);
            %aggiunge gli strati verticali
            interferro_parete = (h-2*s)/(nBarreParete/4 + 1);
            y = s + interferro_parete;
            for i=1:(nBarreParete/4)
%                 obj.addRebar(s+(h-2*s)*(i/(nBarreParete/4)),4*AreaBarra);
                obj.addRebar(y,4*AreaBarra);
                y = y + interferro_parete;
            end
            
            obj.initgeo;
        end
        
        function plotSection(obj)
            %plotta la sezione
            
            for i=1:size(obj.layer,1)
                %plotta i layers
                y = obj.layer(i,1); %ordinata del lembo inferiore del layer
                h = obj.layer(i,2); %altezza del layer
                b = obj.layer(i,3); %larghezza di base del layer
                t = obj.layer(i,4); %larghezza superiore del layer
                Xp = [-b/2 b/2 t/2 -t/2];
                Yp = [y y y+h y+h];
                patch(Xp,Yp,[0.8 0.8 0.8]);
                hold on;
            end
            
            for i=1:size(obj.rebar,1)
                %plotta le barre
                y = obj.rebar(i,1);
                As = obj.rebar(i,2);
%                 line([-As/2*1000 As/2*1000],[y y],'Color','red');
                plot([-As/2*1000 As/2*1000],[y y],'r','LineWidth',3);
            end
			
			hold off;
%             text(0,obj.H/2,num2str(size(obj.rebar,1)));
            axis equal;
        end
        
        
    end
    
end