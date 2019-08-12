function [outputArg1,outputArg2] = CPP_readBumpyJamoutput(inputArg1,inputArg2)
%%
% clear 
% close all
% clc

path = ['/Users/philipwang/Dropbox/Yale/C++/Bumpy/output_6_16_1_1_1.xy'];

fileID = fopen(path,'r');
formatSpec = '%f';
A = textscan(fileID,formatSpec,7,'Delimiter','\n');
Nc = A{1}(1);
n = A{1}(2);
N = A{1}(3);
mu = A{1}(4);
Lx = A{1}(5);
Ly = A{1}(6);
gam = A{1}(7);

bumpy_format = "";
type = "%f, ";
for i=1:n
%     if i == n
%         type = "%f ";
%     end
    bumpy_format = strcat(bumpy_format,type);
end

formatSpec = '%f, %f, %f, %f, ';
formatSpec = formatSpec + bumpy_format + bumpy_format + "%d";

B = textscan(fileID,formatSpec,Nc,'Delimiter','\n');

Dn = B{1}(:)';
x = B{2}(:)';
y = B{3}(:)';
th = B{4}(:)';
r_shape = [];
th_shape = [];
for i=5:5+n-1
    r_shape = [r_shape,B{i}(:)];
end
for i=5+n:5+2*n-1
    th_shape = [th_shape,double(B{i}(:))];
end
Cn = B{end}(:);

plot_bumpy(Nc,n,N,Dn,x,y,th,r_shape,th_shape,Lx,Ly)

end

