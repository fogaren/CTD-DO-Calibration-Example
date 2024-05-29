
% Files with oxygen values output in voltages with default hysteresis
% correction 

a1 = 4; a2 = 5; a3 = 6; % Time lags to try 

files = ls('*_hyst.cnv'); % Depends on file naming convention 
cast = [];
for i = 1:height(files)
    cast = readSBScnv(files(i,:));
    if max(cast.pm) > 250 % Only look at casts deeper than 250 dbar 
        align_CTD_Year9(cast,a1,a2,a3);
    end
end

%%
function align_CTD_Year9(cast,a1,a2,a3) 

SP = gsw_SP_from_C(cast.c0mScm,cast.t090C,cast.pm);
SA = gsw_SA_from_SP(SP,cast.pm,cast.lon,cast.lat);
CT = gsw_CT_from_t(SA,cast.t090C,cast.pm);
cast.pt = gsw_pt_from_CT(SA,CT); 

CastString = ['Cast: ' cast.source(9:end-9)]; % Depends on file naming convention 

cast.DO_a1v = SBE_alignCTDW(cast.sbeox0V, a1, 1/24 );
cast.DO_a2v = SBE_alignCTDW(cast.sbeox0V, a2, 1/24 );
cast.DO_a3v = SBE_alignCTDW(cast.sbeox0V, a3, 1/24 );

pres = 25; % To remove surface noise from plots 

figure('Position',[ 100 50 1200 600]);
subplot(1,6,1)
plot(cast.pt(cast.pm >= pres),cast.pm(cast.pm >= pres))
axis ij
ylabel('Pressure')
xlabel('pot. temp (\circC)')
set(gca,'Fontsize',11)
title('Temp')

subplot(1,6,2)
plot(cast.sbeox0V(cast.pm >= pres),cast.pm(cast.pm >= pres))
axis ij
xlabel('Oxygen (V)')
ylabel('Pressure')
set(gca,'Fontsize',11)
title('Oxygen')

subplot(1,6,3)
plot(cast.sbeox0V(cast.pm >= pres),cast.pt(cast.pm >= pres))
ylabel('Temp')
xlabel('DO (V)')
title('0 sec align')
set(gca,'Fontsize',11)

subplot(1,6,4)
plot(cast.DO_a1v(cast.pm >= pres),cast.pt(cast.pm >= pres))
ylabel('Temp')
xlabel('DO (V)')
title([num2str(a1) ' sec align'])
set(gca,'Fontsize',11)

subplot(1,6,5)
plot(cast.DO_a2v(cast.pm >= pres),cast.pt(cast.pm >= pres))
ylabel('Temp')
xlabel('DO (V)')
title([num2str(a2) ' sec align'])
set(gca,'Fontsize',11)

subplot(1,6,6)
plot(cast.DO_a3v(cast.pm >= pres),cast.pt(cast.pm >= pres))
ylabel('Temp')
xlabel('DO (V)')
title([num2str(a3) ' sec align'])
sgtitle(CastString)
set(gca,'Fontsize',11)

end