%
% Written by Adam Schonewille, University of Toronto, Canada.
% April 2020 (C)

clear all
clc

[file_name file_path]=uigetfile({'*.txt','Data File (.TXT)'},'Select Data','multiselect','off');
% file_name=sort(file_name);
[file_name2 file_path2]=uiputfile('*.txt','Save as formatted Datafile',file_path);
% lps=questdlg('How many loops?','Loops','Forever','None','Other','Forever');
% switch lps
%     case 'Forever'
%         loops=65535;
%     case 'None'
%         loops=1;
%     case 'Other'
%         loops=inputdlg('Enter number of loops? (must be an integer between 1-65535)        .','Loops');
%         loops=str2num(loops{1});
% end
% 
% delay=inputdlg('What is the delay time? (in seconds)        .','Delay');
% delay=str2num(delay{1});
% dly=questdlg('Different delay for the first image?','Delay','Yes','No','No');
% if strcmp(dly,'Yes')
%     delay1=inputdlg('What is the delay time for the first image? (in seconds)        .','Delay');
%     delay1=str2num(delay1{1});
% else
%     delay1=delay;
% end
% dly=questdlg('Different delay for the last image?','Delay','Yes','No','No');
% if strcmp(dly,'Yes')
%     delay2=inputdlg('What is the delay time for the last image? (in seconds)        .','Delay');
%     delay2=str2num(delay2{1});
% else
%     delay2=delay;
% end

numVariables = 4;

% h = waitbar(0,['0% done'],'name','Progress') ;
fid=fopen(strcat(file_path,file_name),'r');
fid2=fopen(strcat(file_path2,file_name2),'w');
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    values = tline(~isspace(tline));
    noSpacesAt = (~isspace(tline));
%     startIndex = 1;
%     count = 0;
    rEdge = 0;
    fEdge = 0;
    tempValuesLengthVector = [];
    for n = 1:length(noSpacesAt)-1
%         count = count+1;
        if ((noSpacesAt(n) - noSpacesAt(n+1)) > 0)
            % falling edge
            fEdge = n;
            dist = fEdge - rEdge;
            tempValuesLengthVector = [tempValuesLengthVector, dist]; 
        end
        if ((noSpacesAt(n) - noSpacesAt(n+1)) < 0)
            % Rising edge
            rEdge = n;
            if ( length(tempValuesLengthVector) == (numVariables-1) )
                dist = length(noSpacesAt) - rEdge;
                tempValuesLengthVector = [tempValuesLengthVector, dist];
            end
        end 
    end
    x = values(1:tempValuesLengthVector(1));
    y = values(tempValuesLengthVector(1)+1:tempValuesLengthVector(1)+tempValuesLengthVector(2));
    z = values(tempValuesLengthVector(1)+tempValuesLengthVector(2)+1:tempValuesLengthVector(1)+tempValuesLengthVector(2)+tempValuesLengthVector(3));
    mag = values(tempValuesLengthVector(1)+tempValuesLengthVector(2)+tempValuesLengthVector(3)+1:end);    
    Data = [str2num(x),str2num(y),str2num(z),str2num(mag)];
    fprintf(fid2,'%.18f\t%.18f\t%.18f\t%.18f\n',Data);
    disp(tline)
    
end
fclose(fid);
fclose(fid2);
% waitbar(i/length(file_name),h,[num2str(round(100*i/length(file_name))),'% done']) ;

return;

for i=1:length(file_name)
    if strcmpi('gif',file_name{i}(end-2:end))
        [M  c_map]=imread([file_path,file_name{i}]);
    else
        a=imread([file_path,file_name{i}]);
        [M  c_map]= rgb2ind(a,256);
    end
    if i==1
        imwrite(M,c_map,[file_path2,file_name2],'gif','LoopCount',loops,'DelayTime',delay1)
    elseif i==length(file_name)
        imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay2)
    else
        imwrite(M,c_map,[file_path2,file_name2],'gif','WriteMode','append','DelayTime',delay)
    end
    waitbar(i/length(file_name),h,[num2str(round(100*i/length(file_name))),'% done']) ;
end

close(h);
msgbox('Finished Successfully!')