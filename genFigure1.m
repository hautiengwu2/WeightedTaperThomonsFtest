close all ; clear

N = 12199 ; %20199 ;
NW = 433 ; %433 ;
W = NW ./ N ;

NoDPSWF = 2*NW ;


t1 = tic ;
if NoDPSWF > 2*NW
    [v, lambda] = dpss(N, NW, NoDPSWF) ; % get all DPSS
else
    [v, lambda] = dpss(N, NW) ; % get all DPSS
end
toc(t1) 


 

%% prepare some normalizations/parameters

UP = 1 ;

c = zeros(NoDPSWF, 1) ;
for jj = 1: NoDPSWF
    c(jj) = sum(v(:, jj)) ;
end

c = c./norm(c) ;


if mod(N*UP, 2) == 0 % even
    xi = (-N*UP/2: N*UP/2-1)' / (N*UP) ;
else % odd (symmetric wrt 0)
    xi = (-(N*UP-1)/2 : (N*UP-1)/2)' / (N*UP) ;
end

SHIFT = exp(sqrt(-1)*pi*xi*(N-1)) ;
MID = find(xi == 0) ;
[~, WIDX] = min(abs(xi-W)) ;


% this is the f coordinate associated with Bk; that is, the final
% oscillation point
%wk = MID + round(acos(Bk)/2/pi /(xi(2)-xi(1))) ;



 

%% evaluate desired quantities
% prepare space to save
DPSWF = zeros(N*UP, NoDPSWF) ;
dDPSWF = zeros(N*UP, NoDPSWF) ;
V = zeros(N*UP, 1) ; % using DPSWF, the original Thomson's
RHO = zeros(N*UP, 1) ;


% record hinge points
Ukat0 = zeros(NoDPSWF, 1) ;
UkatEND = zeros(NoDPSWF, 1) ;
dUkat0 = zeros(NoDPSWF, 1) ;
dUkatEND = zeros(NoDPSWF, 1) ;
Uksup = zeros(NoDPSWF, 1) ;


for jj = 1: NoDPSWF 

    const = (-1).^(jj-1) ;

    if mod(jj-1, 2) == 1
        const = const * 1i ;
    end

    vhat = fftshift(fft(v(:,jj), N*UP)) ;
    DPSWF(:, jj) = const * vhat .* SHIFT ; % ignore e^{i2\pi f(N-1)/2} factor. The same
    dDPSWF(:, jj) = derivative(xi, DPSWF(:, jj)) ; 

    if jj <= NW*2 - ceil((5/(pi^2))*log(100*NW+25)*log(N))
        V = V + c(jj) * DPSWF(:, jj) ;
    end

    Ukat0(jj) = DPSWF(MID, jj) ;
    UkatEND(jj) = DPSWF(end, jj) ;
    dUkat0(jj) = dDPSWF(MID, jj) ;
    dUkatEND(jj) = dDPSWF(end, jj) ;
    Uksup(jj) = max(abs(DPSWF(:, jj))) ;

end
RHO = RHO / (NW*2 - ceil((5/(pi^2))*log(100*NW+25)*log(N))) ;

fprintf(['V(0)=' num2str(V(MID)),'\n']) ;



%% generate oversampled version over [-1/2, 1/2]
UP = 19 ;

if mod(N*UP, 2) == 0 % even
    xi2 = (-N*UP/2: N*UP/2-1)' / (N*UP) ;
else % odd (symmetric wrt 0)
    xi2 = (-(N*UP-1)/2 : (N*UP-1)/2)' / (N*UP) ;
end

MID2 = find(xi2 == 0) ;
SHIFT2 = exp(1i*pi*xi2*(N-1)) ;

DPSWF2 = zeros(N*UP, size(v, 2)) ;
dDPSWF2 = zeros(N*UP, size(v, 2)) ;
V2 = zeros(N*UP, 1) ;
RHO2 = zeros(N*UP, 1) ;


for jj = 1: 2*NW 
    const = (-1).^(jj-1) ;

    if mod(jj-1, 2) == 1
        const = const * sqrt(-1) ;
    end

    vhat2 = fftshift(fft(v(:,jj), N*UP)) ;
    DPSWF2(:, jj) = const * vhat2 .* SHIFT2 ;
    tmp = derivative(xi2, DPSWF2(:, jj)) ; 
    dDPSWF2(:, jj) = tmp ;

    if jj <= NW*2 - ceil((5/(pi^2))*log(100*NW+25)*log(N))
        V2 = V2 + c(jj) * DPSWF2(:, jj) ;
        RHO2 = RHO2 + abs(DPSWF2(:, jj)).^2 ;
    end

end
RHO2 = RHO2 / (NW*2 - ceil((5/(pi^2))*log(100*NW+25)*log(N))) ;

fprintf(['V(0)=' num2str(V2(MID2)),'\n']) ;



L2NORM = sqrt((xi2(2)-xi2(1)) * sum(abs(V2).^2)) ;
L2NORMU = sqrt((xi(2)-xi(1)) * sum(abs(V).^2)) ;
fprintf(['L2 norm of V and V2 = ' num2str(L2NORM) ' ' num2str(L2NORMU) '\n']) ;





% prepare bound for small f
UPPER = nan(size(xi2)) ;
for kk = 1: floor(NW/log(NW))
    idx = find(xi2 > kk/N & xi2 <= (kk+1)/N) ;
    UPPER(idx) = sqrt(N)/pi/kk ;
end





%%
scrsz = get(0,'ScreenSize') ;

h1 = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)*2.5/4]) ;
tiledlayout(2,1,'TileSpacing','tight','Padding','tight')

nexttile 
plot(xi2, abs(V2), 'color', [.7 .7 .7], 'linewidth', 1) ; hold on;
plot(xi, abs(V), 'k', 'linewidth', 2) ;  
plot([-W -W], [min(abs(V)) max(abs(V))], 'b', 'linewidth', 2) ;
plot([W W], [min(abs(V)) max(abs(V))], 'b', 'linewidth', 2) ;
SIG = 1./abs(sin(pi*xi)) ; SIG(SIG>110) = nan ;
plot(xi, SIG, 'r--', 'linewidth', 2) ;
xlim([-3*W 3*W]) ; ylabel('$|\hat{w}|$', 'interpreter', 'latex') ; set(gca,'xtick', []) ;
set(gca,'fontsize', 18) ; ylim([-3 sqrt(N)+1])

nexttile
plot(xi2, log10(abs(V2)), 'color', [.7 .7 .7], 'linewidth', 1) ; hold on ;
plot(xi2, log10(UPPER), 'r', 'linewidth', 3) ; hold on;
plot(xi, log10(abs(V)), 'k', 'linewidth', 3) ; hold on ;
plot([-W -W], [-17 max(log10(abs(V)))], 'b', 'linewidth', 2) ;
plot([W W], [-17 max(log10(abs(V)))], 'b', 'linewidth', 2) ;
axis tight ; xlim([-3*W 3*W]) ; ylabel('$\log_{10}(|\hat{w}|)$', 'interpreter', 'latex') ; 
xlabel('frequency ($f$)', 'interpreter', 'latex') ;
set(gca,'fontsize', 18) ; 

exportgraphics(h1, 'VplotRough.pdf', 'ContentType', 'vector', 'BackgroundColor','none');


h1 = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)*2.5/4]) ;
tiledlayout(2,1,'TileSpacing','tight','Padding','tight')

nexttile 
plot(xi2, abs(V2), 'color', [.7 .7 .7], 'linewidth', 1) ; hold on;
plot(xi, abs(V), 'k', 'linewidth', 2) ;  
xlim([-0.06*W 0.06*W]) ; ylabel('$|\hat{w}|$', 'interpreter', 'latex') ; set(gca,'xtick', []) ;
set(gca,'fontsize', 18) ; ylim([-5 sqrt(N)+5])

nexttile
plot(xi2, log10(abs(V2)), 'color', [.7 .7 .7], 'linewidth', 1) ; hold on ;
plot(xi, log10(abs(V)), 'k', 'linewidth', 3) ;
axis tight ; xlim([-0.06*W 0.06*W]) ; ylabel('$\log_{10}(|\hat{w}|)$', 'interpreter', 'latex') ; 
set(gca,'fontsize', 18) ;  
xlabel('frequency ($f$)', 'interpreter', 'latex') ;

exportgraphics(h1, 'VplotRough2.pdf', 'ContentType', 'vector', 'BackgroundColor','none');


 