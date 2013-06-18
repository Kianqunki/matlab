clear all;
close all;

dir_list=ls;
temp = size(dir_list);
len = temp(1);
line_size = temp(2);
clear temp;

for i=3:len-2
   entry = dir_list(i,:);
   loc = strfind(entry,'jpg');
   if loc
      loc_m = strfind(entry,'m');
      loc_s = strfind(entry,'s');
      if loc_s-loc_m ~= 10
          temp=strcat(entry(1:loc_m+2),'.000000',entry(loc_s:loc+2));
          movefile(entry,temp);
          disp(strcat(entry,' -> ', temp));
          clear temp;
      end
      
   end
    
    
    
    
end

