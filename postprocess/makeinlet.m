str0 = './';
str1 = 'Re_45_';
str2 = '_u.txt';
str3 = '_v.txt';
str4 = '_inlet.dat';
for i = 1:9
    stri = int2str(i);
    A = importdata(strcat(str0,str1,stri,str2));
    YY = A(:,1);
    VXX = A(:,2);
    B = importdata(strcat(str0,str1,stri,str3));
    VYY = B(:,2);
    VX = [flipud(VXX);VXX];
    VY = [-flipud(VYY);VYY];
    Y = [-flipud(YY);YY];
    output = [Y,VX,VY];
    filename = strcat(str1,stri,str4);
    fileID = fopen(filename,'w');
    fprintf(fileID,'%5i\n',length(Y));
    fprintf(fileID,'%23.16f %23.16f %23.16f\n',output');
    fclose(fileID);
end

A = importdata('./Re_46_0_u.txt');
YY = A(:,1);
VXX = A(:,2);
B = importdata('./Re_46_0_v.txt');
VYY = B(:,2);
VX = [flipud(VXX);VXX];
VY = [-flipud(VYY);VYY];
Y = [-flipud(YY);YY];
output = [Y,VX,VY];
fileID = fopen('Re_46_0_inlet.dat','w');
fprintf(fileID,'%5i\n',length(Y));
fprintf(fileID,'%23.16f %23.16f %23.16f\n',output');
fclose(fileID);
