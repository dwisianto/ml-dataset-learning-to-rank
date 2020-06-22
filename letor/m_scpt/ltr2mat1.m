%
% For a given rank text file, conver the file into MPM readable format
%   - 0 qid:3 1:56 2:204 3:204 4:142 5:140 6:189 7:182 8:183 
% Note:
%   - missing value is NULL
%   rankTextToMpmMatFile('/media/sdb1/db/rankPreference/', 'test_txt_aggregate.txt', 'test_txt_aggregate','/media/sdb1/gb/po/pancancer_po/s0901/Do/');
%
%
%
function [dat lbl datMat lblMat] = ltr2mat1(in_loc,in_name,out_loc,out_name)                                                                        


% Read data   %%%%%%%%%%%%%%%%%%%%%%%%%
% Open file & get the first line
disp('- Read in text file');
filename=[in_loc,in_name];
fid=fopen(filename);
line=fgetl(fid);
fclose(fid);

tmp=textscan(line,'%s');
n=numel(tmp{1}); % the length of tmp gives how many columns (fields)

%%
% Now read the file
fid=fopen(filename,'r');
% First column is gene name, rest are RNAseq values
fmt=repmat('%s\t', [1,n-1]);  
fmt=['%f\t', fmt];
dat=textscan(fid,fmt,'EmptyValue',NaN);
fclose(fid);


%% 1st col = rank label
l=dat{1};

%%
% 2nd col = question ID   FORMAT: qid:xx    
% Note: use CStr2String instead of sprintf if file is large
str=sprintf('%s;', dat{2}{:});  % convert to a single string
qid=strread(str, 'qid:%f;');    % get numeric part only


%%
% number of feature
iCtr=1;
for i=3:numel(dat)  % for each data columns, get the data (2nd numeric value after ':')
    if( isempty(strfind(dat{i}{1},'#')) == 0 ) 
        fea_len = iCtr;
        break;
    end
    iCtr=iCtr+1;
end



%%
cnt=1;
x = zeros(length(qid), fea_len);
d = zeros(length(qid), fea_len);
for i=3:fea_len  % for each data columns, get the data (2nd numeric value after ':')
    str=sprintf('%s;', dat{i}{:});  % convert to a single string
    [x(:,cnt) d(:,cnt)] = strread(str,'%f:%f;',length(qid)+1);  % get the 2nd numeric part only 
    cnt=cnt+1;
end 

clear dat; clear str; clear x;


%%%%%%%%%%%%%%%%%% Convert it to MPM format %%%%%%%%%%%%%%%%%%%%%%
disp('- Convert to MPM format');
qid_uniq     = unique(qid);
qid_uniq_len = length(qid_uniq);

% allocate cell 
% Group data & label based on each qid
dat={};
lbl={};
for i=1:qid_uniq_len
    idx=find(qid==qid_uniq(i));
    dat{i}=d(idx,:);
    lbl{i}=l(idx,:);
end
clear d; clear l; 


%find thetas using labeled training queries
%thetas = find_thetas_supervised(lbl, dat)';

datMat = []; %zeros(length(qid), fea_len);
lblMat = []; %zeros(length(qid), fea_len);
for iCtr=1:qid_uniq_len
    %iCurrLen = size(dat{iCtr},1);
    datMat= [datMat;dat{iCtr}];
    lblMat = [lblMat;lbl{iCtr}];
end
clear qid;

%%
disp([out_loc, out_name]);
save([out_loc, out_name,'_dat.mat'], 'dat');
save([out_loc, out_name,'_lbl.mat'], 'lbl');
save([out_loc, out_name,'_dat_mat.mat'], 'datMat');
save([out_loc, out_name,'_lbl_mat.mat'], 'lblMat');


 
disp('done');
end