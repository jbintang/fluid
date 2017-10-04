function [] = HUD(t,wt,dt,tp,iter,err,tstring)

io = 1:iter.o;
figure(1)
semilogy(io,err.p,io,err.u,io,err.v)
title(['Galat, ',tstring])
legend('tekanan','kecepatan u','kecepatan v')
xlabel('Iterasi')
ylabel('Galat')
drawnow expose

clc
fprintf('Waktu Tidak Berdimensi : %4.3e\n',t);
fprintf('Total Langkah : %4.3e\n',t/dt);
fprintf('Wall Time (s) : %4.3e\n',wt);
fprintf('Laju Rata-rata (langkah/s) : %4.3e\n\n',t/dt/wt);
fprintf('Galat\n')
fprintf(' p:%4.3e\n',err.p(iter.o));
fprintf(' u:%4.3e\n',err.u(iter.o));
fprintf(' v:%4.3e\n\n',err.v(iter.o));
fprintf('Total Waktu untuk Iterasi Tekanan : %4.3e\n',tp);