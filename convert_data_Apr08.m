load apr08_01


% THIS Bit loads the data for first part of April 08:
% 
% fid=fopen('../Data_EURUSD_2008Q1/eurjpy_20080401.csv');
% EURJPY=textscan(fid,'%*s%s%d%*s%f%f%*d%*d','delimiter',',','treatAsEmpty','NULL');
% % Get time in fractions of days:
% EURJPY{1}=datenum(EURJPY{1},'yyyy-mm-dd HH:MM:SS');
% 
% % Add milliseconds:
% for q=1:length(EURJPY{1})
%    EURJPY{1}(q)=EURJPY{1}(q)+double(EURJPY{2}(q)/(1000*24*60*60));
% end
% 
% fid=fopen('../Data_EURUSD_2008Q1/eurusd_20080401.csv');
% EURUSD=textscan(fid,'%*s%s%d%*s%f%f%*d%*d','delimiter',',','treatAsEmpty','NULL');
% % Get time in fractions of days:
% EURUSD{1}=datenum(EURUSD{1},'yyyy-mm-dd HH:MM:SS');
% 
% % Add milliseconds:
% for q=1:length(EURUSD{1})
%    EURUSD{1}(q)=EURUSD{1}(q)+double(EURUSD{2}(q)/(1000*24*60*60));
% end
% 
% fid=fopen('../Data_EURUSD_2008Q1/usdjpy_20080401.csv');
% USDJPY=textscan(fid,'%*s%s%d%*s%f%f%*d%*d','delimiter',',','treatAsEmpty','NULL');
% % Get time in fractions of days:
% USDJPY{1}=datenum(USDJPY{1},'yyyy-mm-dd HH:MM:SS');
% 
% % Add milliseconds:
% for q=1:length(USDJPY{1})
%    USDJPY{1}(q)=USDJPY{1}(q)+double(USDJPY{2}(q)/(1000*24*60*60));
% end


% THIS Bit loads the data for first part of Jan 08:% 
% 
% fid=fopen('../Data_EURUSD_2008Q1/BNKA-EURUSD-2008-01-1.csv');
% EURUSD=textscan(fid,'%*.s%s%s%*.s%*.s%n%n','delimiter',',','treatAsEmpty','NULL');
% % Get time in fractions of days:
% EURUSD{1}=datenum(EURUSD{1},'yyyymmdd');
% EURUSD{2}=datenum(EURUSD{2},'HH:MM:SS:FFF');
% % Add milliseconds:
% for q=1:length(EURUSD{1})
% EURUSD{1}(q)=EURUSD{1}(q)+EURUSD{2}(q);
% end
% 
% fid=fopen('../Data_EURUSD_2008Q1/BNKA-USDJPY-2008-01-1.csv');
% USDJPY=textscan(fid,'%*.s%s%s%*.s%*.s%n%n','delimiter',',','treatAsEmpty','NULL');
% % Get time in fractions of days:
% USDJPY{1}=datenum(USDJPY{1},'yyyymmdd');
% USDJPY{2}=datenum(USDJPY{2},'HH:MM:SS:FFF');
% % Add milliseconds:
% for q=1:length(USDJPY{1})
% USDJPY{1}(q)=USDJPY{1}(q)+USDJPY{2}(q);
% end
% 
% 



% plot(USDJPY{1},(USDJPY{3}-mean(USDJPY{3}))./std(USDJPY{3}),'r')
% hold on
% plot(EURUSD{1},(EURUSD{3}-mean(EURUSD{3}))./std(EURUSD{3}),'g')
% plot(EURJPY{1},(EURJPY{3}-mean(EURJPY{3}))./std(EURJPY{3}),'c')

% Resample onto uniform time grid:
min_time=min([EURUSD{1}' EURJPY{1}' USDJPY{1}'])

max_time=max([EURUSD{1}' EURJPY{1}' USDJPY{1}'])

reg_grid=min_time:1/(60*10):max_time;

EURUSD_uniform=convert_to_uniform(EURUSD{1},EURUSD{3},reg_grid);
USDJPY_uniform=convert_to_uniform(USDJPY{1},USDJPY{3},reg_grid);
EURJPY_uniform=convert_to_uniform(EURJPY{1},EURJPY{3},reg_grid);

% 
% diff_t=find(diff(EURUSD{1})>0);
% 
% EURUSD_uniform=interp1(EURUSD{1}(diff_t),EURUSD{3}(diff_t),reg_grid,'nearest','extrap');
% 
% diff_t=find(diff(USDJPY{1})>0);
% 
% USDJPY_uniform=interp1(USDJPY{1}(diff_t),USDJPY{3}(diff_t),reg_grid,'nearest','extrap');
% 
% diff_t=find(diff(EURJPY{1})>0);
% 
% EURJPY_uniform=interp1(EURJPY{1}(diff_t),EURJPY{3}(diff_t),reg_grid,'nearest','extrap');

cross_EURJPY=EURUSD_uniform.*USDJPY_uniform;

% plot(reg_grid,(USDJPY_uniform-mean(USDJPY_uniform(~isnan(USDJPY_uniform))))./std(USDJPY_uniform(~isnan(USDJPY_uniform))),'r')
% hold on
% plot(reg_grid,(EURUSD_uniform-mean(EURUSD_uniform))./std(EURUSD_uniform),'g')
% plot(reg_grid,(EURJPY_uniform-mean(EURJPY_uniform))./std(EURJPY_uniform),'c')
% plot(reg_grid,(cross_EURJPY-mean(cross_EURJPY))./std(cross_EURJPY),'c')

reg_grid=reg_grid-min(reg_grid);

%fid=fopen('Data_EURUSD_2008Q1_Shadow_Detail_notitles.CSV');
%EURUSD=textscan(fid,'%*.s%s%n%s%n%n%d%d','delimiter',',','treatAsEmpty','NULL');
% Get time in fractions of days:
%EURUSD{1}=datenum(EURUSD{1},'dd/mm/yyyy HH:MM:SS');
% Add milliseconds:
%for q=1:length(EURUSD{1})
%EURUSD{1}(q)=24*EURUSD{1}(q)+EURUSD{2}(q)/(24*60*60*1000);
%end
%min_t=min(EURJPY{1})
%plot(EURJPY{1}-min_t,EURJPY{3},EURJPY{1}-min_t,EURJPY{4})
%hold on
%min_t=min(EURUSD{1})
%plot(EURUSD{1}-min_t,EURUSD{4},EURUSD{1}-min_t,EURUSD{5})
