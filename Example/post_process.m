% post_processing program for FEMT
% 
% input: PROGRAM.CTR, PROGRAM.OUT
% 
% Last Edited by  GAO HEXUAN at 2016-2-25

%% initial
clear all;
close all;
clc;


%% open files
fprintf('open files\n');

% tmpline = input('Input the program''s name:','s');
% ctr = fopen(strcat(tmpline,'.ctr'), 'r');
out = fopen("../Result/nl_demo/Disp_ans_GALM.out", 'r');
tmpline = fgetl(out);

while ~strncmp(tmpline, ' COORDINATES',12)
    tmpline = fgetl(out);
    if strncmp(tmpline, 'CURRENT_LOAD', 12)
        [~, load_cur] = strtok(tmpline);
        load_cur = str2num(load_cur);
    end
end

%% read coordinates
fprintf('read coordinates\n');
while ~strncmp(tmpline, ' ELEMENTS',9)
    if strncmp(tmpline, ' COORDINATES',12)
        [~,Nnodes] = strtok(tmpline);
        tmpline = fgetl(out);
        for i = 1 : str2num(Nnodes)
            P(i, :) = str2num(fgetl(out));
        end
    end
    tmpline = fgetl(out);
end

[~,Nelements] = strtok(tmpline);

%% read elements
fprintf('read elements\n');
while ~strncmp(deblank(tmpline), '    ELEMENT_MATERIAL',20)
    tmpline = fgetl(out);
    if strncmp(tmpline, '    ELEMENT_NODES',17)
        for i = 1 : str2num(Nelements)
            T(i, :) = str2num(fgetl(out));
        end
    end
end


%% show point
% fprintf('read analysis point\n');
% while ~strncmp(tmpline,'END SOLUTION',12)
%     tmpline = fgetl(ctr);
%     if strncmp(strtrim(tmpline),"S11",3)
%         [~,tmpline] = strtok(tmpline);
%         S11show = str2num(tmpline);
%     elseif strncmp(strtrim(tmpline),"S12",3)
%         [~,tmpline] = strtok(tmpline);
%         S12show = str2num(tmpline);
%     elseif strncmpi(strtrim(tmpline),"temp",3)
%         [~,tmpline] = strtok(tmpline);
%         tempshow = str2num(tmpline);
%     end
% end
% fclose(ctr);

%% read displacement
fprintf('read displacement\n');
while ~strcmp(tmpline, '*END')
    tmpline = fgetl(out);
    if strcmp(tmpline, ' *** NODE_VALUE ***')
            tmpline = fgetl(out);
        for i = 1 : str2num(Nnodes)
            disp(i, :) = str2num(fgetl(out));
        end
        break
    end
end
frewind(out);

%% read stress
% fprintf('read stress\n');
% while tmpline ~= -1
%     tmpline = fgetl(out);
%     if strcmp(tmpline, ' *** NODAL STRESS ***')
%             tmpline = fgetl(out);
%         for i = 1 : str2num(Nnodes)
%             stress(i, :) = str2num(fgetl(out));
%         end
%     end
% end
% fclose(out);

%% check mesh

%   viewmsh(P(:,2:3),T(:,2:end));
 

%% draw displacement
fprintf('calculate displacement\n');
% scale_factor = max(max(abs(P(:,2:end))))/max(max(abs(disp(:,2:end))))*0.1;
% post_P(:,2:3) = scale_factor*disp(:,2:3) + P(:,2:3);
% post_P(:,1) = P(:,1);
% 
% U1 = disp(:,2);
% U2 = disp(:,3);
% Temper = disp(:,4);
% U  = sqrt(disp(:,2).^2+disp(:,3).^2); 
% r = sqrt(P(:,2).^2 + P(:,3).^2);
% n(:,1:2) = [P(:,2)./r,P(:,3)./r];
% Un = n(:,1).*disp(:,2)+n(:,2).*disp(:,3);
% t(:,1:2) = [-P(:,3)./r,P(:,2)./r];
% Ut = t(:,1).*disp(:,2)+t(:,2).*disp(:,3);
U = disp(:,2);
V = disp(:,3);
% Umag = sqrt(U.^2 + V.^2);

x = P(:, 2);
y = P(:, 3);
x1 = x + U;
y1 = y + V;
% x1= post_P(:,2);
% y1= post_P(:,3);
tt = T(:,2:end);
ttp = T(:,2:5);

[m, u] = size(tt);
if u == 6
    tt = [tt(:, 1), tt(:, 4), tt(:, 2), tt(:, 5), tt(:, 3), tt(:, 6)];
elseif u == 8
    tt = [tt(:, 1), tt(:, 5), tt(:, 2), tt(:, 6), tt(:, 3), tt(:, 7), tt(:, 4), tt(:, 8)];
elseif u == 9
    tt = [tt(:, 1), tt(:, 5), tt(:, 2), tt(:, 6), tt(:, 9), tt(:, 8), tt(:, 8), tt(:, 9), ...
          tt(:, 6), tt(:, 3), tt(:, 7), tt(:, 4)];
end

tt = tt';


fprintf('draw displacement\n');

figure
% subplot(2,2,[2,4])
% patch(x(tt), y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;
% patch(-x(tt), y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;
% patch(x(tt), -y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;
% patch(-x(tt), -y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;

% patch(x1(tt), y1(tt), U(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
patch(x1(tt), y1(tt),U(tt));
colormap(jet(50))
axis image;
box on;
colorbar;
% patch(-x1(tt), y1(tt), U(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
% patch(x1(tt), -y1(tt), U(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
% patch(-x1(tt), -y1(tt), U(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
% drawc(x,y,ttp,Pressure);
if load_cur > 0
    title(sprintf('U-disp(Current Load: %.2f)', load_cur))
else
    title('U-disp(Current Load: Nonequilibrium)')
end
% legend('before','after','Location','NorthEast')
% axis image;
% colorbar;

% subplot(2,2,1)
% drawc(x,y,tt,U);
% title('U1');
% 
% subplot(2,2,3)
% drawc(x,y,tt,V);
% title('U2');

% subplot(2,2,3)
% drawc(x,y,tt,Umag);
% title('U');

%% draw temp
% figure(2)
% drawc(x,y,tt,Tmp);

%% draw stress
% fprintf('calculate stress\n');
% S11 = stress(:,2);
% S22 = stress(:,3);
% S12 = stress(:,4);
% S1  = stress(:,5);
% S2  = stress(:,6);

% for i = 1:length(P(:,2))
%     SN = [n(i,1),n(i,2);t(i,1),t(i,2)]*[stress(i,2),stress(i,4);stress(i,4),stress(i,3)]*[n(i,1),n(i,2);t(i,1),t(i,2)]';
%     SNN(i) = SN(1,1);
%     STT(i) = SN(2,2);
%     SNT(i) = SN(1,2);
% end

% Mises  = sqrt(3/2*((S1-(S1+S2)/3).^2+(S2-(S1+S2)/3).^2));

% 
% fprintf('draw stress\n');
% 
% figure(3)
% subplot(2,2,1)
% drawc(x,y,tt,S11,S11show);
% title('S11');
% 
% subplot(2,2,2)
% drawc(x,y,tt,S22);
% title('S22');
% 
% subplot(2,2,[3,4])
% drawc(x,y,tt,S12,S12show);
% title('S12');


% subplot(2,2,4)
% patch(x(tt), y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;
% patch(-x(tt), y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;
% patch(x(tt), -y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;
% patch(-x(tt), -y(tt), 'w','FaceAlpha',.1,'LineWidth',.1);hold on;

% patch(x1(tt), y1(tt), Mises(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
% patch(-x1(tt), y1(tt), Mises(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
% patch(x1(tt), -y1(tt), Mises(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
% patch(-x1(tt), -y1(tt), Mises(tt),'FaceAlpha',.6,'LineWidth',.1);hold on;
% title('Mises');
% colorbar;
% axis image;
% box off;
% [Umax,UI] = max(abs(Mises));
% text(x1(UI),y1(UI),strcat('\leftarrowmax = ',num2str(sign(U(UI))*Umax)));

