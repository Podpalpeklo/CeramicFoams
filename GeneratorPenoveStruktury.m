classdef GeneratorPenoveStruktury < handle
%% GENETRATORPENOVESTRUKTURY Tøída zajišující práci s jádry pro Voroného teselaci ve 3D
%     Verze: 2016-10-23
%     Kontakt: Roman Papšík <roman.papsik@vutbr.cz>

	properties
		velikostBunky = 0
		rozmery = [0 0 0]
		pravidelnost = 0
        pocetBunek = 0
		jadra = []
	end
	
	methods (Access = public)

		function gps = GeneratorPenoveStruktury
		%% GENERATORPENOVESTRUKTURY Konstruktor tøídy (nepotøebuje žádné parametry)
		%     Použití: instance = GeneratorPenoveStruktury

		end

		function VygenerovatJadra(gps, rozmery, velikostBunky, pravidelnost)
		%% VYGENEROVATJADRA Vygenetruje v prostoru jádra dle zadaných parametrù
		%     instance.VYGENEROVATJADRA(rozmery, velikostBunky, pravidelnost)
		%         rozmery - vektor 1x3
		%         velikostBunky - prùmìr koule buòce vepsané
		%         pravidelnost - hodnota v procentech

			gps.rozmery = rozmery;
			gps.velikostBunky = velikostBunky;
			gps.pravidelnost = pravidelnost/100;

			% pod 60% pravidelnosti použít jinou metodu
			if (gps.pravidelnost < 0.6)
				gps.VygenerovatNahodne;
			else
				gps.VygenerovatPravidelne;
				gps.Rozmistit;
			end

			disp('Generování jader:');
			disp(['  * objem [', num2str(gps.rozmery), ']']);
			disp(['  * pravidelnost [', num2str(gps.pravidelnost*100), '%]']);
			disp(['  * velikost buòky [', num2str(gps.velikostBunky), ']']);
		end

		function NacistJadra(gps, varargin)
		%% NACISTJADRA Naète jádra ze zvoleného souboru
		%     instance.NACISTJADRA - zobrazi dialogové okno pro výbìr souboru
		%     instance.NACISTJADRA('C:\soubor.xyz') - naète soubor pøímo

			if nargin == 2
				adresa = fullfile(varargin{1});
			% pokud nebyl zadany soubor, otevre se dialog
			else
				[jmeno_souboru, adresar_souboru] = uigetfile({'*.xyz', 'Souøadnice jader'}, 'Naèíst soubor s jádry');
				% sobor nebyl zvolen
				if isequal(jmeno_souboru, 0) || isequal(adresar_souboru, 0)
					warning('Soubor nebyl vybrán');
					return;
				% pouzije se zvoleny soubor
				else
					adresa = fullfile(adresar_souboru, jmeno_souboru);
				end
			end
			
			gps.rozmery = csvread(adresa, 0, 0, [0 0 0 2]);
			gps.jadra = csvread(adresa, 1, 0);
			gps.jadra = gps.jadra(:, 1:3);
			disp('Naèítání jader:');
			disp(['  * soubor [', adresa, ']']);
		end

		function UlozitJadra(gps, varargin)
		%% ULOZITJADRA Uloží jádra ze zvoleného souboru
		%     instance.ULOZITJADRA - zobrazi dialogové okno pro výbìr souboru
		%     instance.ULOZITJADRA('C:\soubor.xyz') - uloží soubor pøímo

			% oveøení, zda nìjaká jádra existují
			if isempty(gps.jadra)
				warning('Neexistují žádné jádra');
				return;
			end

			% pokud byl zadany soubor, pouzije se
			if nargin == 2
				adresa = fullfile(varargin{1});
			% pokud nebyl zadany soubor, otevre se dialog
			else
				[jmeno_souboru, adresar_souboru] = uiputfile({'*.xyz', 'Souøadnice jader'}, 'Uložit soubor s jádry');
				% sobor nebyl zvolen
				if isequal(jmeno_souboru, 0) || isequal(adresar_souboru, 0)
					warning('Soubor nebyl vybrán');
					return;
				% pouzije se zvoleny soubor
				else
					adresa = fullfile(adresar_souboru, jmeno_souboru);
				end
			end

			dlmwrite(adresa, [gps.rozmery; gps.jadra], 'precision', '%.18g', 'newline', 'pc');
			disp('Ukládání jader:');
			disp(['  * soubor [', adresa, ']']);
		end

		function VytvoritApdl(gps, typ, varargin)
		%% VYTVORITAPDL Odešle jádra aplikaci pro vygenerování kódu APDL
		%     instance.VYTVORITAPDL(vystupniSoubor) - Parametrem je adresa a název výstupního souboru

			% vstupni soubor se vygeneruje automaticky v temp slozce
			vstup = [tempname, '.xyz'];
			dlmwrite(vstup, [(1:length(gps.jadra))', gps.jadra], 'delimiter', ' ', 'precision', '%-.18g', 'newline', 'pc');

			disp('Vytváøení voroného teselace a APDL kódu...');

			if nargin == 3
				vystup = fullfile(varargin{1});
			else
				[jmeno_souboru, adresar_souboru] = uiputfile({'*.inp', 'Kód APDL'}, 'Uložit vstupní soubor pro Ansys');
				% sobor nebyl zvolen
				if isequal(jmeno_souboru, 0) || isequal(adresar_souboru, 0)
					warning('Soubor nebyl vybrán');
					return;
				% pouzije se zvoleny soubor
				else
					vystup = fullfile(adresar_souboru, jmeno_souboru);
				end
			end

			[status, odpoved] = system(['VoroUMTBM.exe ', typ, ' ', num2str(gps.rozmery, 18), ' ', vstup, ' ', vystup]);

			if(status ~= 0)
				warning('VoroUMTBM.exe selhal');
				return;
			end

			disp('APDL kód vytvoøen');
			disp(['  * výstupní soubor [', vystup, ']']);
		end
	end

	methods (Access = protected)

		function VygenerovatPravidelne(gps)
		%% VYGENEROVATPRAVIDELNE Vytvoøí zcela pravidelnou strukturu BCC møížky

			gps.jadra = []; % vyprazdneni existujicich jader
			x_max = gps.rozmery(1);
			y_max = gps.rozmery(2);
			z_max = gps.rozmery(3);

			pocet_x = ceil(x_max / gps.velikostBunky);
			pocet_y = ceil(y_max / gps.velikostBunky);
			pocet_z = ceil(z_max / gps.velikostBunky);

			if (pocet_x*pocet_y*pocet_z) == 0
				error('Pøíliš velká buòka');
			end

			kompenzace_x = (x_max - pocet_x * gps.velikostBunky)/2;
			kompenzace_y = (y_max - pocet_y * gps.velikostBunky)/2;
			kompenzace_z = (z_max - pocet_z * gps.velikostBunky)/2;

			for j = 0:pocet_y
				for i = 0:pocet_x
					for k = 0:pocet_z
						gps.jadra(end+1,:) = [gps.velikostBunky*i, gps.velikostBunky*j, gps.velikostBunky*k];
						gps.jadra(end,:) = gps.jadra(end,:) + [kompenzace_x, kompenzace_y, kompenzace_z];
					end
				end
			end

			for j = 1:pocet_y
				for i = 1:pocet_x
					for k = 1:pocet_z
						gps.jadra(end+1,:) = [gps.velikostBunky*i, gps.velikostBunky*j, gps.velikostBunky*k];
						gps.jadra(end,:) = gps.jadra(end,:) - [gps.velikostBunky/2, gps.velikostBunky/2, gps.velikostBunky/2];
						gps.jadra(end,:) = gps.jadra(end,:) + [kompenzace_x, kompenzace_y, kompenzace_z];
					end
				end
			end
		end

		function VygenerovatNahodne(gps)
		%% VYGENEROVATNAHODNE Umisuje náhodnì do prostoru jádra, pokud se v jejich okolí nevyskytují jiné
			gps.jadra = []; % pro jistotu vyprazdnit
			x_max = gps.rozmery(1);
			y_max = gps.rozmery(2);
			z_max = gps.rozmery(3);

			pocet_x = ceil(x_max / gps.velikostBunky);
			pocet_y = ceil(y_max / gps.velikostBunky);
			pocet_z = ceil(z_max / gps.velikostBunky);

			pocet_jader = ((pocet_x+1)*(pocet_y+1)*(pocet_z+1))+(pocet_x*pocet_y*pocet_z);

			i = 0;
			while i < pocet_jader % vytvaret body, do dosazeni pozadovaneho poctu
				x = rand * gps.rozmery(1);
				y = rand * gps.rozmery(2);
				z = rand * gps.rozmery(3);
				if(gps.OvereniVzdalenosti([x, y, z])) % prijmout pokud nema jine jadra ve svem okoli
					gps.jadra(i+1,:) = [x, y, z];
					i = i + 1;
				end
			end
		end

		function Rozmistit(gps)
		%% ROZMISTIT Vychýlí jádra z pozic urèených zcela pravidelnou strukturou
		
			for i=1:length(gps.jadra)
				gps.jadra(i,1) = gps.jadra(i,1) + gps.velikostBunky*(1-gps.pravidelnost)*rand*sin(rand*2*pi)*cos(rand*2*pi); %nahodny posun jadra v ramci kulove plochy
				gps.jadra(i,2) = gps.jadra(i,2) + gps.velikostBunky*(1-gps.pravidelnost)*rand*cos(rand*2*pi)*cos(rand*2*pi); %nahodny posun jadra v ramci kulove plochy
				gps.jadra(i,3) = gps.jadra(i,3) + gps.velikostBunky*(1-gps.pravidelnost)*rand*sin(rand*2*pi);                %nahodny posun jadra v ramci kulove plochy
			end
			for i=1:length(gps.jadra)
				for j=1:length(gps.rozmery)
					if gps.jadra(i,j) > gps.rozmery(j)
						gps.jadra(i,j) = gps.rozmery(j);
					end
					if gps.jadra(i,j) < 0
						gps.jadra(i,j) = 0;
					end
				end
			end
		end

		function d_min = MinimalniVzdalenost(gps)
		%% MINIMALNIVZDALENOST Vypoèítá nejmenší vzdálenost mezi jádry pøi dokonalém uspoøádání
			pocet_x = ceil(gps.rozmery(1) / gps.velikostBunky);
			pocet_y = ceil(gps.rozmery(2) / gps.velikostBunky);
			pocet_z = ceil(gps.rozmery(3) / gps.velikostBunky);
			pocet_jader = pocet_x*pocet_y*pocet_z;

			 % Nejmenší dosažitelná vzdálenost mezi jádry v prostoru
			d_min = nthroot(3*sqrt(3)/4/pocet_jader, 3);
		end

		function prijmout = OvereniVzdalenosti(gps, kandidat)
		%% OVERENIVZDALENOSTI Zjistí, zda v daném okolí bodu nejsou jiné jádra
		%     instance.OVERENIVZDALENOSTI([x,y,z]) - parametrem je vektor 3 souøadnic v prostoru
			if isempty(gps.jadra)
				prijmout = true;
				return;
			end
			cx = kandidat(1);
			cy = kandidat(2);
			cz = kandidat(3);
			blizke = find( sqrt(power(gps.jadra(:,1)-cx,2) + power(gps.jadra(:,2)-cy,2) + power(gps.jadra(:,3)-cz,2)) < gps.MinimalniVzdalenost, 1);
			
			if isempty(blizke)
				prijmout = true;
			else
				prijmout = false;
			end
		end

	end
end