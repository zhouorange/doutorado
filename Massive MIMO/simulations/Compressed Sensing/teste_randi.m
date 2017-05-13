 for idx=1:1:10000
a=randi(2048,1,16);
if(~isempty(find(a==2048))) 
    fprintf(1,'Achei!!!\n');
    find(a==2048)
    break 
end
 end