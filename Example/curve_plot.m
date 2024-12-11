%% initial
clear all;
close all;
clc;

%% open files
fprintf('open files\n');

out = fopen("../Result/nl_demo/Load_Disp_curve.txt", 'r');
tmpline = fgetl(out);

id = 1;
while(tmpline ~= -1)
    if(length(tmpline) > 0)
        tmpnum = str2num(tmpline);
        x(id) = tmpnum(1);
        y(id) = tmpnum(2);
        id = id + 1;
    end
    tmpline = fgetl(out);
end
start_label = 16;
len = length(x);
figure
% plot(x(1,start_label:len), y(1,start_label:len), Color='r');hold on;
plot(x(1,1:len), y(1,1:len), Color='r');hold on;
plot(x(1,start_label:len), y(1,start_label:len), LineStyle="none", Marker="o", MarkerEdgeColor='g');hold on;
plot(x(1,1:start_label-1), y(1,1:start_label-1), LineStyle="none", Marker="o",MarkerEdgeColor='b');hold on;
title("屈曲平衡路径");
ylabel("端部载荷p/Pa");
xlabel("最大水平位移U_m/m", "Interpreter","tex");
fclose(out);