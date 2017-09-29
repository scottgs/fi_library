function [ data ] = Grant_loadcsvfile( filename )

%CYGWIN
%  sed -i '1d' filename 
%  sed -e s/[a-zA-Z0-9/_.]*.tif,//g -i filename
%or all
%  sed -i '1d' *
%  sed -e s/[a-zA-Z0-9/_.]*.tif,//g -i *

    data = csvread(filename);

%     fid = fopen( filename );
%     
%     n = 0;
%     tline = fgetl(fid);
%     while ischar(tline)
%       tline = fgetl(fid);
%       n = n+1;
%     end
%         
%     fclose(fid);
%     
%     fid = fopen( filename );
%     
%     data = zeros(n-1,23);
%     
%     l=fgetl(fid);
%     ctr = 1;
%     while 1
%         s = fgets(fid);
%         if(isempty(s)==1)
%            break; 
%         end
%         [~,s] = strtok(s,',');
%         [class,s] = strtok(s,',');
%         [conf,s] = strtok(s,',');
%         [label,s] = strtok(s,',');
%         d = [ str2num(class) str2num(conf) str2num(label) ];
%         for i=1:20
%             [v,s] = strtok(s,',');
%             d = [ d str2num(v) ];
%         end
%         data(ctr,:) = d;
%         ctr = ctr + 1;
%     end
%     
%     fclose(fid);
    
end