% Membersihkan Command Window
clc 
% Membersihkan Figure Window
clf 
% Membersihkan Workspace
clear

% Set waktu awal operasi
wt0 = cputime;

%% User Input

% Input jumlah dimensi yang dipakai, jika dua hanya XY
nd = input('Input Jumlah Dimensi (2/3) [3]: '); 
if isempty(nd)  
    % Jika tidak diinput, otomatis akan memakai 3D
    nd = 3;  
end

% Input panjang kotak pada sumbu x
Lx = input('Input Panjang Dimensi Arah X [1]: ');
if isempty(Lx)
    % Jika tidak diinput, otomatis akan bernilai 1
    Lx = 1;
end

% Input panjang kotak pada sumbu y
Ly = input('Input Panjang Dimensi Arah Y [1]: ');
if isempty(Ly)
    % Jika tidak diinput, otomatis akan bernilai 1
    Ly = 1;
end

if nd == 3
    % (Hanya untuk 3D) Input panjang kotak pada sumbu z 
    Lz = input('Input Panjang Dimensi Arah Z [1]: ');
    if isempty(Lz)
        % Jika tidak diinput, otomatis akan bernilai 1
        Lz = 1;
    end
end

% Input penggunaan baffle atau tidak
baf = input('Apakah menggunakan baffle? Y/N [Y]: ','s');
if isempty(baf)
    % Jika tidak diinput, otomatis akan dipakai
    baf = 'Y';
end

% (Hanya untuk pemakaian baffle) Input posisi baffle
if strcmp(baf,'Y') == 1
    % Input posisi baffle dari dinding sebelah kiri (West)
    Lbw = input('Input Jarak Baffle Dari Dinding Kiri (Kelipatan 0.1) [0.2]: ');
    if isempty(Lbw)
        % Jika tidak diinput, otomatis akan bernilai 0.2
        Lbw = 0.2;
    end
    % Input posisi baffle dari dinding sebelah bawah (South untuk 2D)
    Lbb = input('Input Jarak Baffle Dari Dinding Bawah (Kelipatan 0.1) [0.2]: ');
    if isempty(Lbb)
        % Jika tidak diinput, otomatis akan bernilai 0.2
        Lbb = 0.2;
    end
    % Input posisi baffle dari dinding sebelah atas (North untuk 2D)
    Lbt = input('Input Jarak Baffle Dari Dinding Atas (Kelipatan 0.1) [0.2]: ');
    if isempty(Lbt)
        % Jika tidak diinput, otomatis akan bernilai 0.2
        Lbt = 0.2;
    end
    % Input tebal baffle
    hb = input('Input Tebal Baffle [0.01]: ');
    if isempty(hb)
        % Jika tidak diinput, otomatis akan bernilai 0.01
        hb = 0.01;
    end 
end

if strcmp(baf,'Y') == 1
    % (Hanya jika baffle digunakan) Input sumber kecepatan antara Kipas 
    % atau Surface-Driven
    vel = input('Input Sumber Kecepatan SD/F [F]: ','s');
    if isempty(vel)
        % Jika tidak diinput, otomatis akan memakai kipas
        vel = 'F';
        % Karena tidak ada variasi arah kipas, maka nilai variasinya kosong
        Vv = 0;
    end
else
    % (Hanya jika baffle tidak digunakan) Sumber kecepatan otomatis berupa 
    % Surface-Driven
    vel = 'SD';
end

if strcmp(vel,'SD') == 1
    % (Hanya untuk Surface-Driven) Posisi Sumber Kecepatan
    if nd == 3
        % (Hanya untuk 3D) Input posisi SD antara Top, Bottom, North, South, East
        % West Wall
        V = input('Input Posisi Sumber Kecepatan T/B/N/S/E/W [N]: ','s');
        if isempty(V)
            % Jika tidak diinput, otomatis memakai Top Wall
            V = 'T';
        end 
        if strcmp(V,'T') == 1 || strcmp(V,'B') == 1
            % (Hanya untuk SD Top atau Bottom) Input arah kecepatan antara 
            % sumbu x atau y
            Vv = input('Input Arah Sumber Kecepatan X/Y [X]: ','s');
            if isempty(Vv)
                % Jika tidak diinput, otomatis akan memakai sumbu x
                Vv = 'X';
            end
        
        elseif strcmp(V,'N') == 1 || strcmp(V,'S') == 1
            % (Hanya untuk SD North atau South) Input arah kecepatan antara
            % sumbu x atau z
            Vv = input('Input Arah Sumber Kecepatan X/Z [Z]: ','s');
            if isempty(Vv)
                % Jika tidak diinput, otomatis akan memakai sumbu x
                Vv = 'Z';
            end
        
        elseif strcmp(V,'E') == 1 || strcmp(V,'W') == 1
            % (Hanya untuk SD East atau West) Input arah kecepatan antara sumbu
            % y atau z
            Vv = input('Input Arah Sumber Kecepatan Y/Z [Z]: ','s');
            if isempty(Vv)
                % Jika tidak diinput, otomatis akan memakai sumbu z
                Vv = 'Z';
            end
        end    
    elseif nd == 2
        % (Hanya untuk 2D) Input posisi SD antara North, South, East West Wall
        V = input('Input Posisi Sumber Kecepatan N/S/E/W [N]: ','s');
        if isempty(V)
            % Jika tidak diinput, otomatis akan memakai North Wall
            V = 'N';
        end
    end
else
    % (Hanya untuk kipas) Input posisi kipas dari dinding bawah (South
    % untuk 2D)
    V = input('Input Posisi Kipas Dari Dinding Bawah (Kelipatan 0.1) [0.8]: ');
    if isempty(V)
        % Jika tidak diinput, otomatis akan bernilai 0.8
        V = 0.8;
    end
    Vv = 0;
end

% Input arah sumber kecepatan (positif atau negatif)
arahv = input('Input Arah Sumber Kecepatan 1/-1 [1]: ');
if isempty(arahv)
    % Jika tidak diinput, otomatis akan bernilai positif
    arahv = 1;
end

% Input penggunaan fillet atau tidak
fil = input('Input Fillet Y/N [N]: ','s');
if isempty(fil)
    % Jika tidak diinput, otomatis fillet akan dipakai
    fil = 'N';
end

% Input nilai bilangan Reynolds
Re = input('Input Bilangan Reynolds [100]: ');
if isempty(Re)
    % Jika tidak diinput, otomatis akan bernilai 100
    Re = 100;
end

if nd == 3
    % (Hanya untuk 3D) Input penggunaan gravitasi
    g = input('Input Percepatan Gravitasi [-1]: ');
    if isempty(g)
        % Jika tidak diinput, otomatis akan bernilai -1
        g = -1;
    end
end

% Input jumlah volume kontrol pada sumbu x
ncx = input('Input Banyak Volume Kontrol Arah X (Kelipatan 10) [10]: ');
if isempty(ncx)
    % Jika tidak diinput, otomatis akan bernilai 10
    ncx = 10;
end
if strcmp(baf,'Y') == 1
    % (Hanya jika baffle dipakai) Volume kontrol sumbu x ditambah 1 untuk
    % mengakomodasi baffle
    ncx = ncx+1;
end

% Input jumlah volume kontrol pada sumbu y
ncy = input('Input Banyak Volume Kontrol Arah Y (Kelipatan 10) [10]: ');
if isempty(ncy)
    % Jika tidak diinput, otomatis akan bernilai 10
    ncy = 10;
end

if nd == 3
    % (Hanya untuk 3D) Input jumlah volume kontrol pada sumbu z
    ncz = input('Input Banyak Volume Kontrol Arah Z (Kelipatan 10) [10]: ');
    if isempty(ncz)
        % Jika tidak diinput, otomatis akan bernilai 10
        ncz = 10;
    end
end

% Input besar time-step yang digunakan
dt = input('Input Besar Time Step [0.000125]: ');
if isempty(dt)
    % Jika tidak diinput, otomatis akan bernilai 0.000125
    dt = 0.000125;
end

% Input nilai dari SOR
a = input('Input Koefisien SOR [1.1]: ');
if isempty(a)
    % Jika tidak diinput, otomatis akan bernilai 1.1
    a = 1.1;
end

% Input jumlah step yang diperlukan untuk mengupdate HUD
update = input('Input Jumlah Step untuk Update HUD [100]: ');
if isempty(update)
    % Jika tidak diinput, otomatis akan bernilai 100
    update = 100;
end

% Input jumlah step maksimum sebelum program diberhentikan
tmax = input('Input Jumlah Step Maksimum [1000000]: ');
if isempty(tmax)
    % Jika tidak diinput, otomatis akan bernilai 1000000
    tmax = 1000000;
end

% Input jumlah iterasi maksimum untuk kalkulasi persamaan tekanan poisson
ipmax = input('Input Jumlah Iterasi Tekanan Maksimum [1000]: ');
if isempty(ipmax)
    % Jika tidak diinput, otomatis akan bernilai 1000
    ipmax = 1000;
end

%% Inisialisasi
nx = ncx+2; % jumlah kontrol volume yang dipakai sumbu x (bersama ghost cell)
ny = ncy+2; % jumlah kontrol volume yang dipakai sumbu y (bersama ghost cell)

if strcmp(baf,'Y') == 1
    % (Hanya untuk pemakaian baffle)
    hx = Lx/(ncx-1); % lebar kontrol volume sumbu x
    dxbb = hx-(hb/2); % lebar kontrol volume di bagian kanan dan kiri baffle
    dxb = hb; % jarak antara u pada kontrol volume baffle
    ibw = Lbw/Lx*(ncx-1); % grid point u baffle dari batas west
    if nd == 3
        % (Hanya untuk 3D)
        ibb = Lbb/Lz*ncz; % grid point v baffle dari batas bottom 
        ibt = Lbt/Lz*ncz; % grid point v baffle dari batas top
    else
        % (Hanya untuk 2D)
        ibs = Lbb/Ly*ncy; % grid point baffle dari batas south 
        ibn = Lbt/Ly*ncy; % grid point baffle dari batas north
    end
else
    % (Hanya jika tidak memakai baffle)
    hx = Lx/(ncx); % lebar kontrol volume arah x    
end

hy = Ly/(ncy); % lebar kontrol volume sumbu y
if nd == 3
    % (Hanya untuk 3D)
    nz = ncz+2; % jumlah kontrol volume yang dipakai sumbu z (bersama ghost cell)
    hz = Lz/(ncz); % lebar kontrol volume sumbu z
end

if strcmp(vel,'F') == 1
    % (Hanya untuk penggunaan kipas)
    if nd == 3
        % (Hanya untuk 3D)
        ifan = V/Lz*ncz; % grid point baffle dari batas Bottom
    else
        % (Hanya untuk 2D)
        ifan = V/Ly*ncy; % grid point baffle dari batas South
    end
else
    % Jika tidak digunakan kipas, maka tidak ada jarak kipas (untuk
    % mempermudah mengatur kondisi batas)
    ifan = 0; 
end

if nd == 2
    % Pembuatan string yang berisi bilangan Reynolds dan Jumlah Kontrol
    % Volume untuk 2D
    tstring = ['Re = ', num2str(Re),', Jumlah Sel = ',num2str(ncx),'x',num2str(ncy)];
else
    % Pembuatan string yang berisi bilangan Reynolds dan Jumlah Kontrol
    % Volume untuk 3D
    tstring = ['Re = ', num2str(Re),', Jumlah Sel = ',num2str(ncx),'x',num2str(ncy),'x',num2str(ncz)];
end
            
                
%% Inisialisasi Matriks
if nd == 2
   %% Untuk 2D
    % Inisialisasi matriks tekanan
    p = zeros(nx,ny); 
    % Inisialisasi matriks kecepatan u
    u = zeros(nx,ny); 
    % Inisialisasi matriks kecepatan u
    v = zeros(nx,ny); 
    
    % Inisiliasi matriks kecepatan u untuk nilai prediksi
    us = u;
    % Inisiliasi matriks kecepatan v untuk nilai prediksi
    vs = v;
    
    % Inisialisasi matriks tekanan untuk menyimpan hasil iterasi sebelumnya
    p_old = p;
    % Inisialisasi matriks kecepatan u untuk menyimpan hasil iterasi sebelumnya
    u_old = u; 
    % Inisialisasi matriks kecepatan v untuk menyimpan hasil iterasi sebelumnya
    v_old = v; 
    % Inisialisasi matriks untuk nilai RHS pada persamaan tekanan poisson
    Q = zeros(nx-1,ny-1); 
    
    % Inisialisasi matriks untuk jarak antar u pada sumbu x (dxu)
    dx1 = zeros(nx,ny); 
    if strcmp(baf,'Y') == 1
        % (Hanya untuk pemakaian baffle)
        for i = 1:nx
            if i == 1+ibw || i == 1+ibw+2
                % Jarak antar u pada bagian kiri dan kanan baffle
                dx1(i,:) = dxbb;
            elseif i == 1+ibw+1
                % Jarak antar u pada baffle
                dx1(i,:) = dxb;
            else
                % Jarak antar u selain pada baffle dan sekitarnya
                dx1(i,:) = hx;
            end
        end
    else
        % (Tanpa pemakaian baffle) Jarak antar u untuk 1 kotak 
        dx1(:,:) = hx;
    end
    
    % Jarak antar pusat volume kontrol pada sumbu x (dxp)
    dx2 = (dx1(1:end-1,1:ny) + dx1(2:end,1:ny))/2;
    
    % Inisialisasi matriks untuk jarak antar v (dyv)
    dy1 = zeros(nx,ny);
    % Jarak antar v untuk 1 kotak 
    dy1(:,:) = hy;
    
    % Jarak antar pusat volume kontrol pada sumbu y (dyp)
    dy2 = (dy1(1:nx,1:end-1) + dy1(1:nx,2:end))/2;

    % Inisialisasi posisi pusat volume kontrol dengan dinding pada sumbu x
    mx = zeros(1,nx);
    % Inisialisasi posisi pusat volume kontrol dengan dinding pada sumbu y
    my = zeros(1,ny);
    
    % Posisi pertama setelah batas West pada sumbu x
    mx(2) = dx2(1,2)/2;
    % Posisi pertama setelah batas South pada sumbu y
    my(2) = dy2(2,1)/2; 

    for i=3 : nx-1
        % Posisi-posisi antara titik pertama dan terakhir pada sumbu x
        mx(i)=mx(i-1)+dx2(i-1,2);
    end
    for i=3 : ny-1
        % Posisi-posisi antara titik pertama dan terakhir pada sumbu y
        my(i)=my(i-1)+dy2(2,i-1);
    end

    % Posisi terakhir sebelum batas East pada sumbu x
    mx(end) = mx(end-1) + dx2(end,2)/2;
    % Posisi terakhir sebelum batas North pada sumbu y
    my(end) = my(end-1) + dy2(2,end)/2;
    
    % Koefisien Persamaan Tekanan Poisson
    As = zeros(nx,ny);
    As(:,3:ny-1) = dt./(dy2(:,2:ny-2).*dy1(:,3:ny-1));
    An = zeros(nx,ny);
    An(:,2:ny-2) = dt./(dy2(:,2:ny-2).*dy1(:,2:ny-2));
    Aw = zeros(nx,ny);
    Aw(3:nx-1,:) = dt./(dx2(2:nx-2,:).*dx1(3:nx-1,:));
    Ae = zeros(nx,ny);
    Ae(2:nx-2,:) = dt./(dx2(2:nx-2,:).*dx1(2:nx-2,:));
    Ap = -(Aw + Ae + As + An);
    
elseif nd == 3 
   %% Untuk 2D
    % Inisialisasi matriks tekanan
    p = zeros(nx,ny,nz);
    % Inisialisasi matriks kecepatan u
    u = zeros(nx,ny,nz);
    % Inisialisasi matriks kecepatan v
    v = zeros(nx,ny,nz);
    % Inisialisasi matriks kecepatan w
    w = zeros(nx,ny,nz);

    % Inisialisi matriks kecepatan u untuk nilai prediksi
    us = u;
    % Inisialisi matriks kecepatan v untuk nilai prediksi
    vs = v;
    % Inisialisi matriks kecepatan w untuk hasil nilai prediksi
    ws = w;

    % Inisialisasi matriks tekanan untuk menyimpan hasil iterasi sebelumnya
    p_old = p;
    % Inisialisasi matriks kecepatan u untuk menyimpan hasil iterasi sebelumnya
    u_old = u;
    % Inisialisasi matriks kecepatan v untuk menyimpan hasil iterasi sebelumnya
    v_old = v;
    % Inisialisasi matriks kecepatan w untuk menyimpan hasil iterasi sebelumnya
    w_old = w;
    % Inisialisasi matriks untuk nilai RHS pada persamaan tekanan poisson
    Q = zeros(nx-1,ny-1,nz-1);

    % Inisialisasi matriks untuk jarak antar u pada sumbu x (dxu)
    dx1 = zeros(nx,ny,nz);
    if strcmp(baf,'Y') == 1
        % (Hanya untuk pemakaian baffle)
        for i = 1:nx
            if i == 1+ibw || i == 1+ibw+2
                % Jarak antar u pada bagian kiri dan kanan baffle
                dx1(i,:,:) = dxbb;
            elseif i == 1+ibw+1
                % Jarak antar u pada baffle
                dx1(i,:,:) = dxb;
            else
                % Jarak antar u selain pada baffle dan sekitarnya
                dx1(i,:,:) = hx;
            end
        end
    else
        % (Tanpa pemakaian baffle) Jarak antar u untuk 1 kotak
        dx1(:,:,:) = hx;
    end
    
    % Jarak antar pusat volume kontrol (dxp)
    dx2 = (dx1(1:end-1,:,:) + dx1(2:end,:,:))/2;

    % Inisialisasi matriks untuk jarak antar v (dyv)
    dy1 = zeros(nx,ny,nz);
    % Jarak antar v untuk 1 kotak 
    dy1(:,:,:) = hy;
    
    % Jarak antar pusat volume kontrol (dyp)
    dy2 = (dy1(:,1:end-1,:) + dy1(:,2:end,:))/2;

    % Inisialisasi matriks untuk jarak antar w (dzw)
    dz1 = zeros(nx,ny,nz);
    % Jarak antar w untuk 1 kotak 
    dz1(:,:,:) = hz;
    
    % Jarak antar pusat volume kontrol (dwp)
    dz2 = (dz1(:,:,1:end-1) + dz1(:,:,2:end))/2;

    % Inisialisasi posisi pusat volume kontrol dengan dinding pada sumbu x
    mx = zeros(1,nx);
    % Inisialisasi posisi pusat volume kontrol dengan dinding pada sumbu y
    my = zeros(1,ny);
    % Inisialisasi posisi pusat volume kontrol dengan dinding pada sumbu z
    mz = zeros(1,nz);

    % Posisi pertama setelah batas West pada sumbu x
    mx(2) = dx2(1,2,2)/2;
    % Posisi pertama setelah batas South pada sumbu y
    my(2) = dy2(2,1,2)/2;
    % Posisi pertama setelah batas Bottom pada sumbu z
    mz(2) = dz2(2,2,1)/2; 

    for i=3 : nx-1
        % Posisi-posisi antara titik pertama dan terakhir pada sumbu x
        mx(i)=mx(i-1)+dx2(i-1,2,2);
    end
    for i=3 : ny-1
        % Posisi-posisi antara titik pertama dan terakhir pada sumbu y
        my(i)=my(i-1)+dy2(2,i-1,2);
    end
    for i=3 : nz-1
        % Posisi-posisi antara titik pertama dan terakhir pada sumbu z
        mz(i)=mz(i-1)+dz2(2,2,i-1);
    end
    
    % Posisi terakhir sebelum batas East pada sumbu x
    mx(end) = mx(end-1) + dx2(end,2,2)/2;
    % Posisi terakhir sebelum batas North pada sumbu y
    my(end) = my(end-1) + dy2(2,end,2)/2;
    % Posisi terakhir sebelum batas North pada sumbu z
    mz(end) = mz(end-1) + dz2(2,2,end)/2;

    % Matriks untuk koefisien A pada persamaan tekanan poisson
    As = zeros(nx,ny,nz);
    As(:,3:ny-1,:) = dt./(dy2(:,2:ny-2,:).*dy1(:,3:ny-1,:));
    An = zeros(nx,ny,nz);
    An(:,2:ny-2,:) = dt./(dy2(:,2:ny-2,:).*dy1(:,2:ny-2,:));
    Aw = zeros(nx,ny,nz);
    Aw(3:nx-1,:,:) = dt./(dx2(2:nx-2,:,:).*dx1(3:nx-1,:,:));
    Ae = zeros(nx,ny,nz);
    Ae(2:nx-2,:,:) = dt./(dx2(2:nx-2,:,:).*dx1(2:nx-2,:,:));
    Ab = zeros(nx,ny,nz);
    Ab(:,:,3:nz-1) = dt./(dz2(:,:,2:nz-2).*dz1(:,:,3:nz-1));
    At = zeros(nx,ny,nz);
    At(:,:,2:nz-2) = dt./(dz2(:,:,2:nz-2).*dz1(:,:,2:nz-2));
    Ap = -(Aw + Ae + As + An + At + Ab);
end

% Pengaturan awal penghitung jumlah iterasi
iter.o = 1;
% Pengaturan awal array untuk residu tekanan
res.p(1) = 0;
% Pengaturan awal waktu iterasi tekanan total
tp = 0;

%% Solusi
if nd == 2
    %% Kalkulasi 2 Dimensi
    for t=0:dt:tmax
    
        %% Prediksi
        % Pengaturan kondisi batas
        if strcmp(baf,'Y') == 1
            % Dengan baffle
            [u,v] = setBCsB(u,v,V,ibw,ibs,ibn,ny,vel,arahv,ifan,fil);
        else
            % Tanpa baffle
            [u,v] = setBCs(u,v,V,ny,vel,arahv,fil);
        end
        
        % Penghitungan Hx dan Hy
        [Hx,Hy] = HHP(u,v,dx1,dx2,dy1,dy2,nx,ny,Re);
        if iter.o == 1
            % Untuk iterasi pertama
            us(2:nx-1,1:ny-1) = u(2:nx-1,1:ny-1) + dt*Hx(2:nx-1,1:ny-1);
            vs(1:nx-1,2:ny-1) = v(1:nx-1,2:ny-1) + dt*Hy(1:nx-1,2:ny-1);
        else
            % Untuk iterasi kedua dan seterusnya
            [Hx_old,Hy_old] = HHP(u_old,v_old,dx1,dx2,dy1,dy2,nx,ny,Re);
            us(2:nx-1,1:ny-1) = u(2:nx-1,1:ny-1) ...
                + dt/2*(3*Hx(2:nx-1,1:ny-1) - Hx_old(2:nx-1,1:ny-1));
            vs(1:nx-1,2:ny-1) = v(1:nx-1,2:ny-1) ...
                + dt/2*(3*Hy(1:nx-1,2:ny-1) - Hy_old(1:nx-1,2:ny-1));
        end

        %% Penghitungan Korektor
        for i = 2:nx-1
            for j = 2:ny-1
                % Kalkulasi Q 
                Q(i,j) = ((us(i+1,j) - us(i,j))./dx1(i,j) ...
                    + (vs(i,j+1) - vs(i,j))./dy1(i,j));
            end
        end
     
        % Waktu awal lama iterasi tekanan
        tpp = cputime;
        for c = 1:ipmax
            % Kalkulasi Metode ADI dengan Algoritma Thomas
            [p,res.p(c)] = PoissonADI2P(p,nx,ny,Ap,An,As,Ae,Aw,Q,a);
            if res.p(c) < 1e-5
                % Iterasi berhenti ketika residu di bawah 1e-8
                break
            elseif isnan(res.p(c)) == 1
                % Iterasi berhenti jika nilai residu adalah NaN
                disp('SOLUSI DIVERGEN!')
                beep
                return
            end
        end
    
        if c == ipmax
            disp('Peringatan: Tekanan tidak konevergen')
            disp((res.p(c)))
        end
        % Total waktu iterasi tekanan
        tp = tp + (cputime - tpp);
        
        %% Koreksi
        u(2:nx-1,2:ny-1) = us(2:nx-1,2:ny-1) ...
            - dt./dx2(2:nx-1,2:ny-1).*(p(2:nx-1,2:ny-1) - p(1:nx-2,2:ny-1));
        v(2:nx-1,2:ny-1) = vs(2:nx-1,2:ny-1) ...
            - dt./dy2(2:nx-1,2:ny-1).*(p(2:nx-1,2:ny-1) - p(2:nx-1,1:ny-2));
    
        %% Penghitungan galat
        err.p(iter.o) = norm(p - p_old);
        err.u(iter.o) = norm(u - u_old);
        err.v(iter.o) = norm(v - v_old);
        
        % Simulasi berhenti jika galat untuk u dan v di bawah 1e-8
        if err.u(iter.o) < 1e-8 && err.v(iter.o) < 1e-8
            beep
            wt = cputime - wt0;
            break
        end
    
        p = p - p(5,5); % Normalisasi Tekanan
        
        % Pengaturan nilai tekanan dan kecepatan baru sebagai nilai awal pada iterasi
        % berikutnya
        p_old = p;
        u_old = u;
        v_old = v;
    
        %% Penampilan HUD
        if rem(t,dt*update) == 0 && t ~= 0
            wt = cputime - wt0;            
            HUD(t,wt,dt,tp,iter,err,tstring)
        end
        
        
        % Penambahan nilai iterasi berikutnya (berlebih satu)
        iter.o = iter.o + 1;
    end
    
    %% Post-processing
    % Penyimpanan informasi memory yang dipakai
    [user,sys] = memory;
    HUD(t,wt,dt,tp,iter,err,tstring)
    
    % Pengaturan kembali kondisi batas
    if strcmp(baf,'Y')
        [u,v] = setBCsB(u,v,V,ibw,ibs,ibn,ny,vel,arahv,ifan,fil);
    else
        [u,v] = setBCs(u,v,V,ny,vel,arahv,fil);
    end
    
    % Inisialisasi matriks kecepatan pada pusat volume kontrol
    uc = zeros(nx,ny);
    vc = zeros(nx,ny);

    % Penghitungan nilai kecepatan u dan v pada pusat volume kontrol
    uc(2:end-1,2:end-1) = 0.5*(u(2:end-1,2:end-1) + u(3:end,2:end-1));
    vc(2:end-1,2:end-1) = 0.5*(v(2:end-1,2:end-1) + v(2:end-1,3:end));
    
    % Penghitungan nilai resultan kecepatan pada pusat volume kontrol
    uu = sqrt(uc.^2 + vc.^2);
    
    % Pembuatan matriks posisi
    [MX,MY] = meshgrid(mx,my);

    % Penampilan kontur kecepatan
    figure(2)
    contourf(MX,MY,uu',20)
    colorbar
    axis([0,Lx,0,Ly])
    xlabel('X')
    ylabel('Y')
    title(['Kontur Bidang XY, ',tstring]);

    % Penampilan garis arus kecepatan
    figure(3)
    streamslice(MX,MY,uc',vc',2)
    axis([0,Lx,0,Ly])
    xlabel('X')
    ylabel('Y')
    title(['Garis Arus Bidang XY, ',tstring])

    % Penampilan vektor kecepatan
    figure(4)
    quiver(MX,MY,uc',vc')
    axis([0,Lx,0,Ly])
    xlabel('X')
    ylabel('Y')
    title(['Vektor Kecepatan Bidang XY, ',tstring])

    % Pemberian nama file
    if strcmp(vel,'F')
        % Untuk pemakaian kipas
        fname = ['2DRe',num2str(Re),' nc',num2str(ncy),' dt',num2str(dt),' ',vel,num2str(V),' Baf',baf,'fil',fil,'SOR',num2str(a)];
    else
        % Tanpa pemakaian kipas
        fname = ['2DRe',num2str(Re),' nc',num2str(ncy),' dt',num2str(dt),' ',vel,V,' Baf',baf,'fil',fil,'SOR',num2str(a)];
    end
    
    % Perubahan titik menjadi koma
    fname = strrep(fname,'.',',');

    % Pembuatan folder dengan nama fname
    mkdir(fname);
    % Penulisan alamat folder fname
    pdir=[pwd,'\',fname];
    % Pembukaan folder fname
    cd(pdir);
    % Penyimpanan workspace pada folder fname
    save(fname);
    hold on;
    % Penyimpanan figure
    saveas(1,['Error ',tstring],'jpeg');
    saveas(1,['Error ',tstring],'fig');
    saveas(2,['Contour ',tstring],'jpeg');
    saveas(2,['Contour ',tstring],'fig');
    saveas(3,['Streamline ',tstring],'jpeg');
    saveas(3,['Streamline ',tstring],'fig');
    saveas(4,['Vector ',tstring],'jpeg');
    saveas(4,['Vector ',tstring],'fig');
    % Pengembalian posisi folder ke semula
    cd ..
    hold off;
    
    
else
    %% Kalkulasi 3 Dimensi
    for t=0:dt:tmax
    
        %% Prediksi
        % Pengaturan kondisi batas
        if strcmp(baf,'Y')
            % Dengan baffle
            [u,v,w] = setBCsB3D(u,v,w,V,ibw,ibb,ibt,nx,ny,nz,vel,Vv,arahv,ifan,fil);
        else
            % Tanpa baffle
            [u,v,w] = setBCs3D(u,v,w,V,vel,Vv,arahv,fil);
        end
        
        % Penghitungan Hx, Hy, dan Hz
        [Hx,Hy,Hz] = HHP3D(u,v,w,dx1,dx2,dy1,dy2,dz1,dz2,nx,ny,nz,Re,g);
        if iter.o == 1
            % Untuk iterasi pertama
            us(2:nx-1,1:ny-1,1:nz-1) = u(2:nx-1,1:ny-1,1:nz-1) + dt*Hx(2:nx-1,1:ny-1,1:nz-1);
            vs(1:nx-1,2:ny-1,1:nz-1) = v(1:nx-1,2:ny-1,1:nz-1) + dt*Hy(1:nx-1,2:ny-1,1:nz-1);
            ws(1:nx-1,1:ny-1,2:nz-1) = w(1:nx-1,1:ny-1,2:nz-1) + dt*Hz(1:nx-1,1:ny-1,2:nz-1);
        else
            % Untuk iterasi kedua dan setelahnya
            [Hx_old,Hy_old,Hz_old] = HHP3D(u_old,v_old,w_old,dx1,dx2,dy1,dy2,dz1,dz2,nx,ny,nz,Re,g);
            us(2:nx-1,1:ny-1,1:nz-1) = u(2:nx-1,1:ny-1,1:nz-1) ...
                + dt/2.*(3*Hx(2:nx-1,1:ny-1,1:nz-1) - Hx_old(2:nx-1,1:ny-1,1:nz-1));
            vs(1:nx-1,2:ny-1,1:nz-1) = v(1:nx-1,2:ny-1,1:nz-1) ...
                + dt/2.*(3*Hy(1:nx-1,2:ny-1,1:nz-1) - Hy_old(1:nx-1,2:ny-1,1:nz-1));
            ws(1:nx-1,1:ny-1,2:nz-1) = w(1:nx-1,1:ny-1,2:nz-1) ...
                + dt/2.*(3*Hz(1:nx-1,1:ny-1,2:nz-1) - Hz_old(1:nx-1,1:ny-1,2:nz-1));
        end
    
        %% Penghitungan Korektor
        for i = 2:nx-1
            for j = 2:ny-1
                for k = 2:nz-1
                    % Penghitungan Q
                    Q(i,j,k) = ((us(i+1,j,k) - us(i,j,k))./dx1(i,j,k) ...
                        + (vs(i,j+1,k) - vs(i,j,k))./dy1(i,j,k) ...
                        + (ws(i,j,k+1) - ws(i,j,k))./dz1(i,j,k));
                end
            end
        end
     
        % Waktu awal lama iterasi tekanan
        tpp = cputime;
        for c = 1:ipmax
            % Kalkulasi Metode ADI dengan Algoritma Thomas
            [p,res.p(c)] = PoissonADI2P3D(p,nx,ny,nz,Ap,An,As,Ae,Aw,At,Ab,Q,a);
            if res.p(c) < 1e-5
                % Iterasi berhenti ketika 1e-5
                break
            elseif isnan(res.p(c)) == 1
                % Iterasi berhenti jika nilai residu adalah NaN
                disp('SOLUTION DIVERGED')
                beep
                return
            end
        end
    
    
        if c == ipmax
            disp('Peringatan: Tekanan Tidak Konvergen')
            disp((res.p(c)))
        end
        % Total waktu iterasi tekanan
        tp = tp + (cputime - tpp); 
    
        %% Koreksi
        u(2:nx-1,2:ny-1,2:nz-1) = us(2:nx-1,2:ny-1,2:nz-1) ...
            - dt./dx2(1:nx-2,2:ny-1,2:nz-1).*(p(2:nx-1,2:ny-1,2:nz-1) - p(1:nx-2,2:ny-1,2:nz-1));
        v(2:nx-1,2:ny-1,2:nz-1) = vs(2:nx-1,2:ny-1,2:nz-1) ...
            - dt./dy2(2:nx-1,1:ny-2,2:nz-1).*(p(2:nx-1,2:ny-1,2:nz-1) - p(2:nx-1,1:ny-2,2:nz-1));
        w(2:nx-1,2:ny-1,2:nz-1) = ws(2:nx-1,2:ny-1,2:nz-1) ...
            - dt./dz2(2:nx-1,2:ny-1,1:nz-2).*(p(2:nx-1,2:ny-1,2:nz-1) - p(2:nx-1,2:ny-1,1:nz-2));
        
        %% Penghitungan galat
        err.p(iter.o)=0;
        err.u(iter.o)=0;
        err.v(iter.o)=0;
        err.w(iter.o)=0;
    
        for i=1:nz
            erp = norm(p(:,:,i) - p_old(:,:,i));
            err.p(iter.o) = err.p(iter.o) + erp;
            eru = norm(u(:,:,i) - u_old(:,:,i));
            err.u(iter.o) = err.u(iter.o) + eru;
            erv = norm(v(:,:,i) - v_old(:,:,i));
            err.v(iter.o) = err.v(iter.o) + erv;
            erw = norm(w(:,:,i) - w_old(:,:,i));
            err.w(iter.o) = err.w(iter.o) + erw;
        end
    
        % Simulasi berhenti jika galat untuk u, v, dan w di bawah 1e-8
        if err.u(iter.o) < 1e-8 && err.v(iter.o) < 1e-8 && err.w(iter.o) < 1e-8
            beep
            wt = cputime - wt0;
            break
        end

        p = p - p(5,5,5); % Normalisasi tekanan
    
        p_old = p;
        u_old = u;
        v_old = v;
        w_old = w;
    
        %% Penampilan HUD
        if rem(t,dt*update) == 0 && t ~= 0  
            wt = cputime - wt0;
            HUD3D(t,wt,dt,tp,iter,err,tstring)
        end
        iter.o = iter.o + 1;
    end
    % Penyimpanan informasi memory yang dipakai
    [user,sys] = memory;
    HUD3D(t,wt,dt,tp,iter,err,tstring)
    
    % Pengaturan kembali kondisi batas
    if strcmp(baf,'Y')
        [u,v,w] = setBCsB3D(u,v,w,V,ibw,ibb,ibt,nx,ny,nz,vel,Vv,arahv,ifan,fil);
    else
        [u,v,w] = setBCs3D(u,v,w,V,vel,Vv,arahv,fil);
    end
    
    % Penggabungan dimensi kotak menjadi satu array
    L = [Lx Ly Lz];
    
    % Inisialisasi matriks kecepatan pada pusat volume kontrol
    uc = zeros(nx,ny,nz);
    vc = zeros(nx,ny,nz);
    wc = zeros(nx,ny,nz);

    % Penghitungan nilai kecepatan u, v, dan w pada pusat volume kontrol
    uc(2:end-1,2:end-1,2:end-1) = 0.5*(u(2:end-1,2:end-1,2:end-1) + u(3:end,2:end-1,2:end-1));
    vc(2:end-1,2:end-1,2:end-1) = 0.5*(v(2:end-1,2:end-1,2:end-1) + v(2:end-1,3:end,2:end-1));
    wc(2:end-1,2:end-1,2:end-1) = 0.5*(w(2:end-1,2:end-1,2:end-1) + w(2:end-1,2:end-1,3:end));

    % Pemberian nama file
    if strcmp(vel,'F') == 1
        % Dengan kipas
        fname = ['3DRe',num2str(Re),'nc',num2str(ncx),num2str(ncy),num2str(ncz),'dt',num2str(dt),vel,num2str(V),num2str(Vv),'Baf',baf,'g',num2str(g),'fil',fil,'SOR',num2str(a)];
    else
        % Tanpa kipas
        fname = ['3DRe',num2str(Re),'nc',num2str(ncx),num2str(ncy),num2str(ncz),'dt',num2str(dt),vel,V,Vv,'Baf',baf,'g',num2str(g),'fil',fil,'SOR',num2str(a)];
    end
    
    % Perubahan titik menjadi koma
    fname = strrep(fname,'.',',');

    % Pembuatan folder dengan nama fname
    mkdir(fname);
    % Penulisan alamat folder fname
    pdir=[pwd,'\',fname];
    % Pembukaan folder fname
    cd(pdir);
    % Penyimpanan workspace pada folder fname
    save(fname);
    % Penyimpanan grafik galat
    saveas(1,['Error ',tstring],'jpg');
    saveas(1,['Error ',tstring],'fig');
    cd ..
    
    in = 1;
    while in~=0,
        % Input penampilan grafik yang diinginkan
        disp('Pilih gambar yang diinginkan : ');
        disp('(0) Keluar')
        disp('(1) Vektor Kecepatan 3D')
        disp('(2) Vektor Kecepatan Bidang XY')
        disp('(3) Vektor Kecepatan Bidang XZ')
        disp('(4) Vektor Kecepatan Bidang YZ')
        disp('(5) Garis Arus Bidang XY')
        disp('(6) Garis Arus Bidang XZ')
        disp('(7) Garis Arus Bidang YZ')
        disp('(8) Kontur Bidang XY')
        disp('(9) Kontur Bidang XZ')
        disp('(10) Kontur Bidang YZ')
        in = input('');
        if isempty(in)  
            % Jika tidak diinput, otomatis akan bernilai 0
            in = 0;  
        end
        
        if in==1
            % Penampilan dan penyimpanan vektor kecepatan 3D
            vector3D(mx,my,mz,uc,vc,wc,L,tstring,pdir);
        
        elseif in==2
            % Penampilan dan penyimpanan vektor kecepatan bidang XY
            vectorxy(mx,my,mz,uc,vc,L,tstring,pdir);

        elseif in==3
            % Penampilan dan penyimpanan vektor kecepatan bidang XZ
            vectorxz(mx,my,mz,uc,wc,L,tstring,pdir);
        
        elseif in==4
            % Penampilan dan penyimpanan vektor kecepatan bidang YZ        
            vectoryz(mx,my,mz,vc,wc,L,tstring,pdir);
        
        elseif in==5
            % Penampilan dan penyimpanan garis arus bidang XY
            streamxy(mx,my,mz,uc,vc,L,tstring,pdir);
        
        elseif in==6
            % Penampilan dan penyimpanan garis arus bidang XZ
            streamxz(mx,my,mz,uc,wc,L,tstring,pdir);
            
        elseif in==7
            % Penampilan dan penyimpanan garis arus bidang YZ
            streamyz(mx,my,mz,vc,wc,L,tstring,pdir);
            
        elseif in==8
            % Penampilan dan penyimpanan kontur bidang XY
            contourxy(mx,my,mz,uc,vc,L,tstring,pdir);

        elseif in==9
            % Penampilan dan penyimpanan kontur bidang XZ
            contourxz(mx,my,mz,uc,wc,L,tstring,pdir);
            
        elseif in==10
            % Penampilan dan penyimpanan kontur bidang YZ
            contouryz(mx,my,mz,vc,wc,L,tstring,pdir);
        end
    end
end