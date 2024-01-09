% ########################################################################
%                          +-------------------+
%                          |  M-CHI-N STRISCE  |
%                          +-------------------+
%                          Autore: Andrea Marchi
%                              Versione: 1.4
%                                09/01/2024
%
% In fondo al documento sono presenti il ToDo e il ChangeLog.
% Questa classe è basata sul codice "SeReS_strisce" di Andrea Marchi.
% ########################################################################

classdef MchiNstrisce3 < handle
    properties
        hmax % altezza massima del layer
        layer % matrice dove ogni riga è un layer e le colonne sono: y, h, b_bot, b_top
        rebar % matrice dove ogni riga è una barra e le colonne sono: y, As, mat_id
%         materials % struttura che contiene le funzioni sigma(eps) dei materiali (ogni riga è un materiale differente)
        legameCLS % legame sigma(eps) del CLS
        legameAcc % legame sigma(eps) dell'acciaio
        % proprietà aggiuntive:
        layersProp % ogni riga è un layer, le colonne sono: Area, ygl
        % proprietà dei materiali:
        fc % resistenza calcestruzzo
        fy % resistenza acciaio
        epscy, epscu % deformazione di snervamento e ultima del CLS
        epssy, epssu % deformazione di snervamento e ultima dell'acciaio
        % costanti della sezione:
        Area    % Area della Sezione
        yg      % coordinata del baricentro
        Ix      % Momento di inerzia secondo l'asse orizzontale
        H       % altezza totale della sezione
    end
    methods
        function obj = MchiNstrisce3(hmax,CLS,Acc)
            % costruttore della classe
            obj.hmax = hmax;
            obj.legameCLS = CLS;
            obj.legameAcc = Acc;
        end
        
        function addLayer(obj,h,b,t)
            % aggiunge tanti layer alti massimo h_max fino a creare la patch desiderata
            y = 0;
            for i=1:size(obj.layer,1); y = y + obj.layer(i,2); end
            sum = 0;
            while (h - (sum+obj.hmax) > 0)
                bj = b + (t-b)/(y+h)*(y+sum);
                tj = b + (t-b)/(y+h)*(y+sum+obj.hmax);
                if size(obj.layer,1) > 0; obj.layer = [obj.layer; (obj.layer(size(obj.layer,1),1)+obj.layer(size(obj.layer,1),2)) obj.hmax bj tj];
                else; obj.layer = [obj.layer; 0 obj.hmax b tj]; end % aggiunge il primo layer
                sum = sum + obj.hmax;
            end
            bj = b + (t-b)/(h)*(sum);
            if size(obj.layer,1) > 0; obj.layer = [obj.layer; (obj.layer(size(obj.layer,1),1)+obj.layer(size(obj.layer,1),2)) h-sum bj t];
            else; obj.layer = [obj.layer; 0 h b t]; end % aggiunge il primo layer
        end
        
        function addRebar(obj,y,As)
            % aggiunge uno strato di armature ordinarie
            obj.rebar = [obj.rebar; y As];
            
        end
        
        function defineMaterials(obj, fc, fy, eps_cy, eps_cu, eps_sy, eps_su)
            % aggiunge le proprietà dei materiali
            if nargin < 4
                obj.fc = -abs(fc); % negativo: compressione
                obj.fy = fy;
                obj.epscy = -0.002; % negativo: accorciamento
                obj.epscu = -0.0035;
                obj.epssy = 0.002;
                obj.epssu = 0.03; % il massimo a cui sono arrivato è il 6% (le NTC mi pare prescrivano l'1%)
            else
                obj.fc = -abs(fc);
                obj.fy = fy;
                obj.epscy = -abs(eps_cy);
                obj.epscu = -abs(eps_cu);
                obj.epssy = eps_sy;
                obj.epssu = eps_su;
            end
        end
        
        function initgeo(obj)
            % calcola la geometria della sezione
            % si prende in considerazione il fatto che sia composta tutta dello stesso materiale (in realtà dovrei SEMPRE omogeneizzare i layer [per ogni layer un coeff. di omogenizzazione n])
            obj.Area = 0;
            obj.Ix = 0;
            Sx = 0;
            for i=1:size(obj.layer,1)
                % 1  2  3  4 
                % y, h, b, t, mat_id
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
            % calcola N e M impostata la deformazione al lembo inferiore
            % (eps_b) e al lembo superiore (eps_t)
            
            % forze nei layers:
            Flayer = obj.legameCLS((eps_t-eps_b)*obj.layersProp(:,2)/obj.H+eps_b, obj.epscy,obj.epscu,obj.fc) .* obj.layersProp(:,1);
            % forze nelle barre:
            Frebar = obj.legameAcc((eps_t-eps_b)*obj.rebar(:,1)/obj.H+eps_b, obj.epssy,obj.epssu,obj.fy) .* obj.rebar(:,2);
            
            % contributo del calcestruzzo
            N = sum(Flayer);
            M = sum(Flayer .* obj.layersProp(:,2));
            % contributo dell'acciaio
            N = N + sum(Frebar);
            M = M + sum(Frebar .* obj.rebar(:,1));
            
            M = M - N * obj.yg; % aggiunge il momento dovuto alla forza normale
        end
        
        function [N,M] = integrateStrain(obj,eps_inf,eps_sup)
            % ritorna gli sforzi interni alla sezione discretizzando in
            % maniera incrementale in funzione della discrepanza tra
            % tensione approssimata e tensione reale
            Nint = 0;
            Mint = 0;
            for i=1:size(obj.layer,1)
%                 temp = obj; % fa una copia della sezione per il calcolo del layer (così lo può dividere quanto vuole senza ripercussioni) PARE NON FUNZIONARE
                old_layers = obj.layer;
                old_prop   = obj.layersProp;
                [N,M] = controllaLayer(obj,eps_inf,eps_sup,i);
                obj.layer      = old_layers;
                obj.layersProp = old_prop;
                Nint = Nint + N;
                Mint = Mint + M - N*obj.layersProp(i,2);
            end
            N = Nint;
            M = -Mint; % convenzione italiana
        end
        
        function [M,phi,yn] = yieldingPoint(obj, N)
            % calcola il punto di primo snervamento dato uno sforzo normale N
            % DA MIGLIORARE L'ALGORITMO!!!!!!! (vedi la procedura usata per
            % "ultimatePoint" [di MchiNsteel??])
            
            if N < obj.fc*obj.Area; error('la sezione non può resistere a tale compressione'); end
            if N > obj.fy*sum(obj.rebar(:,2)); error('la sezione non può resistere a tale trazione'); end
            
            % fissato lo eps_inf = epssy cambia eps_top fino a trovare N corretto:
            eps_inf = obj.epssy;
            
            % trova un range per fare la bisezione:
            count = 0;
            eps_supP = obj.epssy; % l'ho scelto io (l'N iniziale può essere solamente minore)
            eps_supN = obj.epssy; % valore Negativo
            while obj.setStrain(eps_inf,eps_supN) > N
                eps_supN = eps_supN - 0.001;
                count = count + 1;
                if count > 1000; error('MchiNstrisce non è arrivato a convergenza'); end
            end
            
            % fa la bisezione:
            err = [100 100 100];
            count = 0;
            while ((abs(err(2)) > abs(0.005 * N)) )%|| (abs(eps_supP-eps_supN) > abs(0.0001*eps_supN))) %+0.000001 serve per non creare errori per N=0
                err = [(obj.setStrain(eps_inf,eps_supP) - N) (obj.setStrain(eps_inf,(eps_supP+eps_supN)/2) - N) (obj.setStrain(eps_inf,eps_supN) - N)];
                if err(1)*err(2) < 0
                    % la risultante è tra err(1) e err(2) (eps_supP è corretto)
                    eps_supN = (eps_supP+eps_supN)/2;
                    
                else
                    % la risultante è tra err(2) e err(3) (eps_supN è corretto)
                    eps_supP = (eps_supP+eps_supN)/2;
                end
                count = count + 1;
                if count > 1000; error('MchiNstrisce non è arrivato a convergenza'); end
            end
            
            eps_sup = (eps_supP+eps_supN)/2;
            
            [~,M] = obj.setStrain(eps_inf,eps_sup);
            phi = (eps_inf-eps_sup)/obj.H; % compressione: negativa (la fibra superiore è compressa)
            yn = -eps_sup/(eps_inf-eps_sup) * obj.H;
            M = -M; % convenzione italiana
        end
        
        function [M,phi,yn] = ultimatePoint(obj, N)
            % calcola il momento e la curvatura ultimi
            % il valore di N è considerato negativo di compressione
            
            if N < obj.fc*obj.Area; error('la sezione non può resistere a tale compressione'); end
            if N > obj.fy*sum(obj.rebar(:,2)); error('la sezione non può resistere a tale trazione'); end
            
            % calcola lo sforzo normale interno della sezione nello stato di
            % rottura limite (armature a eps_su e calcestruzzo a eps_cu) e
            % se è comunque minore dello sforzo agente N (quindi se la
            % sezione è più compressa dello sforzo normale agente) vuol dire
            % che l'asse neutro deve "alzarsi" e quindi deve aumentare la
            % epsilon superiore (aumentare perchè prima era negativa:
            % -0.35%) Questo ultimo caso porta ad una rottura lato acciaio
            % (ovvero che è fissa eps_inf e la bisezione si fa su eps_sup)
            if obj.setStrain(obj.epssu,obj.epscu) > N
                % si ha una rottura lato calcestruzzo, quindi è tenuto fisso eps_sup = eps_cu
                rottura_acciaio = 0;
                eps_sup = obj.epscu; % deformazione ultima del calcestruzzo (negativa: di accorciamento)
            else
                % si ha una rottura lato acciaio, quindi è tenuto fisso eps_inf = eps_su
                rottura_acciaio = 1;
                eps_inf = obj.epssu; % deformazione ultima dell'accaio (positiva: di allungamento)
            end
            
%             if rottura_acciaio; warning('rottura lato acciaio'); end
%             if rottura_acciaio; logFile('MchiN::ultimatePoint: rottura lato acciaio'); end
            
            % fa la bisezione:
            err = [100 100 100]; % variabile di errore per la bisezione: [negativo, medio, positivo]
            count = 0;
            if rottura_acciaio == 0
                % rottura del CLS:
                eps_infP = obj.epssu; % deformazione ultima dell'acciaio (non può essere maggiore)
                eps_infN = obj.epscu; % deformazione ultima del calcestruzzo (non può essere minore)
                while ((abs(err(2)) > abs(0.0001 * N + 1e-18)) || (abs(eps_infP-eps_infN) > abs(0.0001*eps_infN))) % +1e-18 serve per non creare errori per N=0)
                    eps_infM = (eps_infP + eps_infN) / 2; % epsilon medio tra i due estremi (risparmio tre calcoli per ogni iterazione al costo di aggiungere una variabile)
                    err = [(obj.setStrain(eps_infN,eps_sup) - N) (obj.setStrain(eps_infM,eps_sup) - N) (obj.setStrain(eps_infP,eps_sup) - N)];
                    if err(2) > 0
                        % la risultante è tra err(1) e err(2) (eps_infN è corretto)
                        eps_infP = eps_infM;
                    else
                        % la risultante è tra err(2) e err(3) (eps_infP è corretto)
                        eps_infN = eps_infM;
                    end
                    count = count + 1;
                    if count > 10000; error('MchiNstrisce non è arrivato a convergenza'); end
                end
                eps_inf = eps_infM;
            else
                % rottura dell'acciaio:
                eps_supP = 0; % potrebbe essere maggiore, potrebbe essere la deformazione di snervamento dell'acciaio
                eps_supN = obj.epscu; % deformazione ultima del calcestruzzo (non può essere minore)
                while ((abs(err(2)) > abs(0.0001 * N + 1e-18)) || (abs(eps_supP-eps_supN) > abs(0.0001*eps_supN))) % +1e-18 serve per non creare errori per N=0)
                    eps_supM = (eps_supP + eps_supN) / 2; % epsilon medio tra i due estremi (risparmio tre calcoli per ogni iterazione al costo di aggiungere una variabile)
                    err = [(obj.setStrain(eps_inf,eps_supN) - N) (obj.setStrain(eps_inf,eps_supM) - N) (obj.setStrain(eps_inf,eps_supP) - N)];
                    if err(2) > 0
                        % la risultante è tra err(1) e err(2) (eps_supN è corretto)
                        eps_supP = eps_supM;
                    else
                        % la risultante è tra err(2) e err(3) (eps_supP è corretto)
                        eps_supN = eps_supM;
                    end
                    count = count + 1;
                    if count > 9000
                        disp('sta per fallire');
                    end
                    if count > 10000; error('MchiNstrisce non è arrivato a convergenza'); end
                end
                eps_sup = eps_supM;
            end
            
            [~,M] = obj.setStrain(eps_inf,eps_sup);
            phi = (eps_inf-eps_sup)/obj.H; % compressione: negativa (la fibra superiore è compressa) IN REALTà è IL CONTRARIO (NEGATIVA)
            yn = -eps_sup/(eps_inf-eps_sup) * obj.H; % yn definito partendo all'"alto"
            M = -M; % convenzione italiana
        end

        
        function eps = deformazione(obj,eps0,chi,y)
            % funzione di epsilon: def(eps0,chi,y) = eps0 + chi*(y-yg)
            eps = eps0 + chi*(y-obj.yg);
        end
        
        function [f] = F(obj,eps0,chi,N0,M0)
            % calcola le funzioni (FUNZIONE AUSILIARIA)
            [N,M]  = obj.Fint(eps0,chi);
            f(1,1) = N - N0;
            f(2,1) = M - M0;
        end
        
%         function [N,M] = Fint(obj,eps0,chi)
%             %restituisce le forze interne della sezione dato un profilo di
%             %deformazione
%             
%             N = sum(obj.legameCLS(obj.deformazione(eps0,chi,(obj.layersProp(:,2)-obj.yg)), obj.epscy,obj.epscu,obj.fc) .* obj.layersProp(:,1)) + ...
%                 sum(obj.legameAcc(obj.deformazione(eps0,chi,(obj.rebar(:,1)-obj.yg)), obj.epssy,obj.epssu,obj.fy) .* obj.rebar(:,2));
%             M = sum(obj.legameCLS(obj.deformazione(eps0,chi,(obj.layersProp(:,2)-obj.yg)), obj.epscy,obj.epscu,obj.fc) .* (obj.layersProp(:,2)-obj.yg) .* obj.layersProp(:,1)) + ...
%                 sum(obj.legameAcc(obj.deformazione(eps0,chi,(obj.rebar(:,1)-obj.yg)), obj.epssy,obj.epssu,obj.fy) .* (obj.rebar(:,1)-obj.yg) .* obj.rebar(:,2));
%         end
        
        function [N,M] = Fint(obj,eps0,chi)
            % restituisce le forze interne della sezione dato un profilo di
            % deformazione (definito da eps0 e chi)
 
            % ----------------------------------------------------------
            % NUOVA FORMULAZIONE SENZA -yg NELLA FUNZIONE "deformazione"
            % ----------------------------------------------------------
            
            N = sum(obj.legameCLS(obj.deformazione(eps0,chi,obj.layersProp(:,2)), obj.epscy,obj.epscu,obj.fc) .* obj.layersProp(:,1)) + ...
                sum(obj.legameAcc(obj.deformazione(eps0,chi,obj.rebar(:,1)), obj.epssy,obj.epssu,obj.fy) .* obj.rebar(:,2));
            M = sum(obj.legameCLS(obj.deformazione(eps0,chi,obj.layersProp(:,2)), obj.epscy,obj.epscu,obj.fc) .* (obj.layersProp(:,2)-obj.yg) .* obj.layersProp(:,1)) + ...
                sum(obj.legameAcc(obj.deformazione(eps0,chi,obj.rebar(:,1)), obj.epssy,obj.epssu,obj.fy) .* (obj.rebar(:,1)-obj.yg) .* obj.rebar(:,2));
        end
        
        function [M,chi] = curvaMchi2(obj,N,chi_lim)
            % calcola la curva M-chi-N
            
            if nargin < 3; chi_max = (obj.epssu-obj.epscu)/obj.H; else; chi_max = chi_lim; end
            
            N_punti = 100; % numero punti della curva
            delta = 1.e-12; % variazione per il calcolo numerico della derivata parziale
            tol = 1.e-8; % tolleranza per la convergenza
            
            chi = linspace(0,chi_max,N_punti);
            M   = zeros(1,N_punti);
            
            for i=2:(N_punti) % parte da 2 perchè così M(0) = 0 (altrimenti non lo è per problemi numerici...)
                err   = 1;
                eps0  = 0;
                count = 0;
                bisez = 0; % controlla se non bisogna ricorrere alla bisezione
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
                    % cerca il range di laboro per la bisezione:
                    eps0_min = obj.epscu; % non può essere minore di questo
                    eps0_max = obj.epscu;
                    Nint = obj.Fint(eps0_max,chi(i)); % ATTENZIONE: se la curvatura è elevata questo potrebbe essere maggiore di N (in tutti gli altri casi DOVREBBE essere MINORE)
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
                                % C'è QUALCOSA CHE NON VA PROPRIO
                                disp('ERRORE nella ricerca del range per la bisezione (per Nint > N)');
                                logFile('MchiN::curvaMchi2::bisezione: ERRORE nella ricerca del range per la bisezione (per Nint > N)');
                                M(i) = 0;
                                return;
                            end
                        end
                    end
%                     eps0_max = eps0_max + 0.0001; %aumenta eps0_max per motivi numerici (così gli estremi dell'intervallo sono più distanti dal punto)
                    % esegue la bisezione:
                    count = 0;
                    Ntot = [1 1 1];
                    tol_bis = 0.00001; % tolleranza per la bisezione
                    while (norm(Ntot) > (tol_bis * abs(N))) || (abs(abs(eps0_max)-abs(eps0_min)) > (tol_bis * abs(eps0_min)))
                        eps0 = (eps0_min + eps0_max) / 2; % valore medio
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
                            eps0_max = eps0_min; % si INVENTA un risultato (tanto un punto solo toppato non influenza molto la ricerca di chi_y tramite minimi quadrati)
                            Ntot = 0; % USCITA FORZATA
                        end
                    end
                end
                [~,M(i)] = obj.Fint(eps0,chi(i));
            end
            
        end
        
        function [M,chi] = curvaMchi3(obj,N,chi_lim)
            % calcola la curva M-chi-N
            %    inizialmente calcola i punti della curva con il metodo di
            %    Newton-Raphson e se fallisce calcola con la bisezione
            %    Con N=0 da problemi numerici, deve essere sempre <0
            
            if nargin < 3; chi_max = (obj.epssu-obj.epscu)/obj.H; else; chi_max = chi_lim; end
            
            N_punti = 500; % numero punti della curva
            delta = 1.e-12; % variazione per il calcolo numerico della derivata parziale
            tol = 1.e-8; % tolleranza per la convergenza (per il metodo di Newton)
            
            chi = linspace(0,chi_max,N_punti);
            M   = zeros(1,N_punti);
            
            for i=2:(N_punti) % parte da 2 perchè così M(0) = 0 (altrimenti non lo è per problemi numerici...)
                err   = 1; % variabile di errore (per il test di convergenza)
                eps0  = 0; % punto iniziale
                count = 0; % contatore che monitora il caso in cui va in loop infinito
                bisez = 0; % controlla se non bisogna ricorrere alla bisezione
                Nint  = 0;
                deltaNint = 100;
                while (err > tol) || (deltaNint > 0.001*abs(N))
                    Nint_old = Nint;
                    [Nint,~] = obj.Fint(eps0,chi(i));
                    deltaNint = abs(Nint-Nint_old);
                    derivata = (obj.Fint(eps0+delta,chi(i))-Nint)/delta;
                    eps0_old = eps0;
                    eps0 = eps0 + (1/derivata) * (N - Nint);
                    if eps0 ~= 0; err = abs((eps0 - eps0_old)/eps0); else; err = abs(eps0 - eps0_old); end
                    count = count + 1;
                    if count > 100 || derivata == 0
%                         warning('MchiN::curvaMchi2: Newton-Raphson non va bene, calcolo con la bisezione...');
                        logFile('MchiN::curvaMchi3: Newton-Raphson non va bene, calcolo con la bisezione...');
                        bisez = 1;
                        break;
                    end
                end
                if bisez
                    % cerca il range di laboro per la bisezione:
                    % fa la ricerca attraverso il "grid search"
                    for j=0:4
                        temp = linspace(-obj.epssu,obj.epssu,(2^j)*10);
                        Nint = zeros(size(temp));
                        for k=1:size(temp,2)
                            Nint(k) = obj.Fint(temp(k),chi(i));
                        end
                        [~,i_min] = min(Nint);
                        [~,i_max] = max(Nint);
                        if (Nint(i_min)<N && Nint(i_max)>N); break; end %ha trovato il range
                        % altrimenti continua finchè non è discretizzato a sufficienza
                    end
                    if (Nint(i_min)>N || Nint(i_max)<N)
                        % se alla fine del for non ha ancora trovato il
                        % range allora non ci sono speranze neanche per la
                        % bisezione
%                         error('MchiN::curvaMchi2::bisezione: ERRORE nella ricerca del range');
                        logFile('MchiN::curvaMchi3::bisezione: ERRORE nella ricerca del range');
                        % provare ad usare il "gradient descend" dal punto
                        % temp(i_min)
                        % IL TEST SULLA SEZIONE CIRCOLARE FALLISCE A QUESTO
                        % PUNTO PER \nu=0.3
                        M(i) = 0;
                        return;
                    end
                    eps0_min = temp(i_min);
                    eps0_max = temp(i_max);
                    
                    % esegue la bisezione:
                    count = 0;
                    Ntot = [1 1 1];
                    tol_bis = 0.00001; % tolleranza per la bisezione
                    while (norm(Ntot) > (tol_bis * abs(N))) || (abs(abs(eps0_max)-abs(eps0_min)) > (tol_bis * abs(eps0_min)))
                        eps0 = (eps0_min + eps0_max) / 2; % valore medio
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
                            logFile('MchiN::curvaMchi3::bisezione: ERRORE convergenza');
                            error('MchiN::curvaMchi3::bisezione: non è arrivato a convergenza');
                        end
                    end
                end
                [~,M(i)] = obj.Fint(eps0,chi(i));
            end
        end
        
        function [My,chiy,Mu,chiu] = findPoints(obj, N, ~)
            % trova i punti di snervamento e ultimo della sezione dalla
            % curva curvatura-momento.
            
            [~,chi_max] = obj.ultimatePoint(N); % METODO PER STABILIZZARE I RISULTATI NUMERICI
            [M,chi] = obj.curvaMchi3(N,chi_max); % calcola la curva "grezza"
            
            % trova il punto ultimo come il punto in cui la curva comincia
            % ad avere la derivata negativa (i materiali non hanno tratti
            % degradanti a tangente negativa):
%             i = find(diff(M) < 0 ,1); % tronca quando comincia a decrescere CERTE VOLTE FALLISCE PER PROBLEMI NUMERICI ALL'INIZIO
            [~,i] = max(M); % tronca quando raggiunge il massimo
            Mu   = M(i);
            chiu = chi(i);
            % elimina la parte di curva dopo il punto ultimo:
            M   = M(1:i);
            chi = chi(1:i);
            
            % trova il punto di snervamento bilinearizzando la curva:
            P0 = [chiu/3 3/4*Mu]; % parametri iniziali di prova
            model = @(P,x) (P(2)/P(1)*x).*heaviside(P(1)-x) + ((Mu-P(2))/(chiu-P(1))*(x-P(1))+P(2)).*heaviside(x-P(1)); % crea il modello da fittare
            lb = [0 0]; % lower bound (non ci possono essere valori negativi)
            ub = [chiu Mu]; % upper bound (il punto di snervamento non può essere maggiore di quello ultimo)
            options = optimoptions('lsqcurvefit','Display','off');
            Pfit = lsqcurvefit(model,P0,chi,M,lb,ub,options); % fitta il modello minimizzando i quadrati
            chiy = Pfit(1);
            My   = Pfit(2);
            
            % controlla i risultati:
            
            % comandi per i test:
            if nargin > 2
                plot(chi,M); hold on;
                plot(linspace(0,chiu,100),model(Pfit,linspace(0,chiu,100)),'k--');
                plot(chiy,My,'bo');
                plot(chiu,Mu,'ro');
            end
            
        end
        
        function [My,chiy,Mu,chiu] = findPoints2(obj, N, ~)
            % trova i punti trovando prima chi_u e poi discretizzando la
            % curva in un numero fisso di punti
            
            [Mu,chiu] = obj.ultimatePoint(N);
            [M,chi] = obj.curvaMchi3(N,chiu); % calcola la curva
            
            % trova il punto di snervamento bilinearizzando la curva:
            P0 = [chiu/3 3/4*Mu]; % parametri iniziali di prova
            model = @(P,x) (P(2)/P(1)*x).*heaviside(P(1)-x) + ((Mu-P(2))/(chiu-P(1))*(x-P(1))+P(2)).*heaviside(x-P(1)); % crea il modello da fittare
            lb = [0 0]; % lower bound (non ci possono essere valori negativi)
            ub = [chiu Mu]; % upper bound (il punto di snervamento non può essere maggiore di quello ultimo)
            options = optimoptions('lsqcurvefit','Display','off');
            Pfit = lsqcurvefit(model,P0,chi,M,lb,ub,options); % fitta il modello minimizzando i quadrati
            chiy = Pfit(1);
            My   = Pfit(2);
            
            % DEBUGGING
            if nargin > 2
                plot(chi,M); hold on;
                plot(chiy,My,'ob');
                plot(chiu,Mu,'or');
            end
            
        end
        
        function out = test(obj,eps0,chi)
            % esegue dei test sulla sezione
            
            N = zeros(size(eps0));
            M = zeros(size(eps0));
            
            for i=1:size(eps0,2)
                [N(i),M(i)] = obj.Fint(eps0(1,i),chi);
            end
            
            out(:,1) = N;
            out(:,2) = M;
            
        end
        
        function initRect(obj,h,b,As1,As2,c)
            % inizializza una sezione rettangolare
            if nargin < 6; c = 0.03; end
            obj.addLayer(h,b,b); % crea la sezione in CLS
            obj.addRebar(c,As1);
            obj.addRebar(h-c,As2);
            
            obj.initgeo;
        end
        
        function initRect2(obj,h,b,rhotot,rhoc_rhotot,db,c)
            % inizializza una sezione rettangolare (con anche l'armatura di parete)
            if nargin < 7; c = 0.03; end
            Ac = b*h; % area calcestruzzo
            Ab = pi/4 * db^2; % area della barra
            np = max(round((h-2*c)/0.3-1), 0); % numero di armature di parete (non ci può essere in interferro maggiore di 30 cm tra esse) [maggiore uguale a zero..]
            n2 = round(Ac*rhotot/Ab * rhoc_rhotot); % numero di armature superiori (compresse)            
            n1 = round(Ac*rhotot/Ab - n2 - 2*np);  % numero di armature inferiori (tese)
            
            obj.addLayer(h,b,b); % crea la sezione in CLS
            obj.addRebar(c,n1*Ab); % crea l'armatura inferiore
            obj.addRebar(h-c,n2*Ab); % crea l'armatura superiore
            % aggiunge le armature di parete:
            if np > 0
                interferro_laterale = (h-2*c)/(np+1);
                y = c + interferro_laterale;
                for i=1:np
                    obj.addRebar(y,2*Ab);
                    y = y + interferro_laterale;
                end
            end
            
            obj.initgeo; % calcola la geometria della sezione
        end
        
        function initCirc(obj,D,db,rhotot,c)
            % inizializza una sezione circolare
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
            % inizializza una sezione rettangolare cava:
            if nargin < 7; c = 0.03; end
            AreaBarra = pi/4 * db^2;
            nBarre = round(rhotot*(8*s*(h+b-2*s))/(pi*db^2));
            interferro = 4*(h+b-2*s)/nBarre;
            nBarreParete = round(4*(h-2*s)/interferro); % round(8/pi * ((h-2*s)*s)/(db^2));
            nBarreStrati = round((nBarre-nBarreParete)/4);
            
            % crea la sezione in CLS:
            obj.addLayer(s,b,b);
            obj.addLayer(h-2*s,2*s,2*s);
            obj.addLayer(s,b,b);
            % aggiunge gli strati orizzontali:
            obj.addRebar(c,AreaBarra*nBarreStrati);
            obj.addRebar(s-c,AreaBarra*nBarreStrati);
            obj.addRebar(h-s+c,AreaBarra*nBarreStrati);
            obj.addRebar(h-c,AreaBarra*nBarreStrati);
            % aggiunge gli strati verticali:
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
            % plotta la sezione
            
			% plotta i layers:
            for i=1:size(obj.layer,1)
                y = obj.layer(i,1); % ordinata del lembo inferiore del layer
                h = obj.layer(i,2); % altezza del layer
                b = obj.layer(i,3); % larghezza di base del layer
                t = obj.layer(i,4); % larghezza superiore del layer
                Xp = [-b/2 b/2 t/2 -t/2];
                Yp = [y y y+h y+h];
                patch(Xp,Yp,[0.8 0.8 0.8]);
                hold on;
            end
            
			% plotta le barre:
            for i=1:size(obj.rebar,1)
                y = obj.rebar(i,1);
                As = obj.rebar(i,2);
                if(0)
                    % stampa le armature come "barrette":
                    line([-As/2*1000 As/2*1000],[y y],'Color','red');
                    plot([-As/2*1000 As/2*1000],[y y],'r','LineWidth',3);
                else
                    % stampa le armature come "cerchietti":
                    r = sqrt(As/pi);
                    theta = linspace(0,2*pi,100);
                    plot(r*cos(theta), y + r*sin(theta), 'r');
                    text(r+0.05, y, ['As=' num2str(As) 'm^2'], 'color','red');
                end
            end
			
			hold off;
%             text(0,obj.H/2,num2str(size(obj.rebar,1)));
            axis equal;
        end
        
        
        function dividiLayer(obj,i)
            % divide in due il layer i (e poi ovviamente aggiorna anche il vettore layerProp
            
            y = obj.layer(i,1);
            h = obj.layer(i,2);
            b = obj.layer(i,3);
            t = obj.layer(i,4);
            
            % modifica il vecchio layer dividendolo nella metà inferiore:
            obj.layer(i,:) = [y, 0.5*h, b, 0.5*(b+t)];
            % aggiunge il nuovo layer nella metà superiore:
            newLayer = [y+0.5*h, 0.5*h, 0.5*(b+t), t];
            obj.layer = [obj.layer(1:i,:); newLayer; obj.layer(i+1:end,:)];
            
            % ricalcola le proprietà per i due layers appena creati:
            Areal = 0.5 * obj.layer(i,2) * (obj.layer(i,3) + obj.layer(i,4));
            ygl   = obj.layer(i,1) + obj.layer(i,2)/3 * (1+obj.layer(i,4)/(obj.layer(i,3)+obj.layer(i,4)));
            obj.layersProp(i,:) = [Areal, ygl];
            Areal = 0.5 * obj.layer(i+1,2) * (obj.layer(i+1,3) + obj.layer(i+1,4));
            ygl   = obj.layer(i+1,1) + obj.layer(i+1,2)/3 * (1+obj.layer(i+1,4)/(obj.layer(i+1,3)+obj.layer(i+1,4)));
            obj.layersProp = [obj.layersProp(1:i,:); [Areal, ygl]; obj.layersProp(i+1:end,:)];
        end
        
        
    end
    
end




function [N,M] = controllaLayer(sez, eps_inf, eps_sup, i)
    % funzione che controlla il layer e se non è verificato lo divide
    % ricorsivamente finchè la distribuzione degli sforzi approssimati
    % all'ordine zero non sono come quelli reali
    
    tol_layers = 0.01; % tolleranza sul controllo dei layers
    
    sigma_inf = sez.legameCLS(eps_inf+(eps_sup-eps_inf)/sez.H*(sez.layer(i,1)), sez.epscy,sez.epscu,sez.fc); %tensione al lembo inferiore del layer
    sigma_sup = sez.legameCLS(eps_inf+(eps_sup-eps_inf)/sez.H*(sez.layer(i,1)+sez.layer(i,2)), sez.epscy,sez.epscu,sez.fc); %tensione al lembo superiore del layer
    sigma_bar = sez.legameCLS(eps_inf+(eps_sup-eps_inf)/sez.H*(sez.layersProp(i,2)), sez.epscy,sez.epscu,sez.fc); %tensione nel baricentro del layer (quella usata nel calcolo approssimato)
    
    % calcola i parametri di controllo dell'errore. dipende dalla tensione e
    % dallo spessore della sezione in quel punto
    err_inf = (sigma_bar - sigma_inf) * sez.layer(i,3);
    err_sup = (sigma_bar - sigma_sup) * sez.layer(i,4);
    
    % controllo sul layer (t e b devono essere diverse da zero [almeno uno dei due])
%     if sigma_bar == 0; err = abs(err_inf)+abs(err_sup); else; err = (abs(err_inf)+abs(err_sup))/(sigma_bar*0.5*(sez.layer(i,3)+sez.layer(i,4))); end %questa condizione non mi piace molto perché è più restrittiva vicino all'asse neutro, dove non ce n'è bisogno...
    if sez.fc == 0; err = abs(err_inf)+abs(err_sup); else; err = (abs(err_inf)+abs(err_sup))/(sez.fc); end %ACCROCCO
    if abs(err) > tol_layers
        % non ha passato il controllo e divide i layers
        sez.dividiLayer(i); %divide il layer in due
        [N1,M1] = controllaLayer(sez,eps_inf,eps_sup,i);
        [N2,M2] = controllaLayer(sez,eps_inf,eps_sup,i+1); % RISUDDIVISI I LAYERS NON è PIù i+1..........!!!!!!
        N = N1 + N2;
        M = M1 + M2 + (N1-N2)*sez.layer(i,2); % l'altezza del layer diviso è hi/2
    else
        N = sigma_bar * sez.layersProp(i,1);
        M = 0;
    end
end




% TODO:
% -----
%   1)  nella funzione "initgeo" va inizializzata la sezione OMOGENEIZZATA
%   2)  modificare "yieldingPoint" 
%   3)  la curvatura in curvaMchi è negativa se comprime le fibre superiori
%   4)  controllare la sezione CIRCOLARE!! (mi pare che veniva diversa da VcaSLU, controllare...)
%   5)  



% ChangeLog:
% ----------
%	13/01/2021: Ho modificato le funzioni "ultimatePoint" eliminando la ricerca iniziale del range e "Fint" eliminando il -yg nel calcolo dellal deformazione
%   16/01/2021: Ho eliminato "findEquilibrium" e "curvaMchi" che non erano usate da tempo
%   17/01/2021: Ho inserito la funzione "curvaMchi3" che migliora "curvaMchi2"
%   21/01/2021: Creato la versione 2 inserendo "integrateStrain", "controllaLayer" e "dividiLayer" 
%   22/02/2021: Aggiunta la visualizzazione delle armature come "cerchietti" nella funzione "plotSection()"
%   09/01/2024: Creato i progetto GitHub e messi insieme gli script che avevo
%   