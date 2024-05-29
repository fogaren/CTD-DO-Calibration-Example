function dcc=import_dcc(filename)

%read_dcc reads downcast or upcast calibrated ctd files for ctd_gui
% modified 9/24/2010 to add secondary O2 sensor

fid=fopen(filename);
lin=fgetl(fid);
n=length(lin);
i=1;
while i<n
    if strcmp(lin(i:i),':')
        dcc.station=str2num(lin(i+1:n));
        i=n;
    else
        i=i+1;
    end
end

lin=fgetl(fid);
if ~isempty(str2num(lin(11:18)))
    dcc.lat=str2num(lin(11:18));
    dcc.lon=str2num(lin(32:40));
end
dcc.date=lin(49:54);
% if length(lin) < 67 
% keyboard
% end
if length(lin) > 54             % added tkm
    dcc.time=lin(63:67);
else
    dcc.time='00:00';
end
lin=fgetl(fid);
space=strfind(lin,' ');
c=0;
for i=1:length(space)-1
    if length(lin(space(i):space(i+1)))>2
        c=c+1;
parms{c}=strtrim(lin(space(i):space(i+1)));
    end
end

parms{c+1}=strtrim(lin(space(end)+1:end));

nparm=length(parms);
form=[];
for i=1:nparm
    form=[form,'%f'];
end

q=fscanf(fid,form,[nparm inf]);
q=q';
fclose(fid);
woce=int2str(q(:,nparm));


for i=1:nparm
   switch parms{i}
   case 'Pres'
       dcc.prs=q(:,i);
       dcc.woce.prs=woce(:,i);
   case 'T90(1)'
       dcc.t901=q(:,i);
       dcc.woce.t901=woce(:,i);
   case 'T90(2)'
       dcc.t902=q(:,i);
       dcc.woce.t902=woce(:,i);
   case 'Sal(1)'
       dcc.sal1=q(:,i);
       dcc.woce.sal1=woce(:,i);
   case 'Sal(2)'
       dcc.sal2=q(:,i);
       dcc.woce.sal2=woce(:,i);
   case 'OxCur'
       dcc.oxcr=q(:,i);
       dcc.woce.oxcr=woce(:,i);
   case 'dovdt'                              %
       dcc.dovdt=q(:,i);
       dcc.woce.dovdt=woce(:,i);
   case 'OXYG'
       dcc.ox=q(:,i);
       dcc.woce.ox=woce(:,i);
   case 'OXYG(ml/l)'
       dcc.ox=q(:,i);
       dcc.woce.ox=woce(:,i);
   case 'uOXY(um/kg)'
       dcc.oxumkg=q(:,i);
       dcc.woce.oxumkg=woce(:,i);
   case 'OxCur2'
       dcc.oxcr2=q(:,i);
       dcc.woce.oxcr2=woce(:,i);
   case 'dovdt2'                              %
       dcc.dovdt2=q(:,i);
       dcc.woce.dovdt2=woce(:,i);
   case 'OXYG2'
       dcc.ox2=q(:,i);
       dcc.woce.ox2=woce(:,i);    
   case 'OXYG2(ml/l)'
       dcc.ox2=q(:,i);
       dcc.woce.ox2=woce(:,i); 
   case 'uOXY2(um/kg)'
       dcc.ox2umkg=q(:,i);
       dcc.woce.ox2umkg=woce(:,i); 
   case 'pH'
       dcc.ph=q(:,i);
       dcc.woce.ph=woce(:,i);
   case 'Trans'
       dcc.tran=q(:,i);
       dcc.woce.tran=woce(:,i);
   case 'Flur'
       dcc.flu=q(:,i);
       dcc.woce.flu=woce(:,i);
   case 'flECO-AFL'
       dcc.flc=q(:,i);                       %5/3/2009 -- renamed flc ( in case there are two
       dcc.woce.flc=woce(:,i);              % as in RU09
   case 'PAR'
       dcc.par=q(:,i);
       dcc.woce.par=woce(:,i);
   case 'BAT'                              %5/3/2009 -- added for RU09
       dcc.bat=q(:,i);
       dcc.woce.bat=woce(:,i);
   case 'SPAR'
       dcc.spar=q(:,i);
       dcc.woce.spar=woce(:,i);
   case 'Turbidity'
       dcc.tur=q(:,i);
       dcc.woce.tur=woce(:,i);
   case 'seaTurbMtr'                         %added by tkm 7/7/2005
       dcc.tur=q(:,i);
       dcc.woce.tur=woce(:,i);
   case 'Altimeter'
       dcc.alt=q(:,i);
       dcc.woce.alt=woce(:,i);
   case 'nscans'
       dcc.ns=q(:,i);
   end
end

