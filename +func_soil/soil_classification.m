%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Input:
%sand: sand fraction [0..1], (nx1)
%clay: Clay fraction[0..1], (nx1)
%T: type of classification: 'usda'
%Varargin{1}: PLOT, BOOLEAN: 1: plots the usda texture triangle plus data
%points
%
%Output:
% SC:    CellArray of soil class strings, (nx1)
% silt:  Silt fraction [%], (nx1)
% SCINT: Array of soil class integers, (nx1)
% SC2:  CellArray of soil class strings in short notation, e.g. S for SAND (nx1)
%
%Example:
%soil_classification(sand,clay, 'usda')
%soil_classification([0.5;0.3;0.1;0],[0.1;0.2;0.4;0.8], 'usda',1)
%--------------------------------------------------------------------------
%disp('-------------------------------------------------------------------')
%disp('                       U S D A soil Classification                 ')
%disp('-------------------------------------------------------------------')
function [SC, silt, SCINT, SC2]=soil_classification(sand, clay, T,varargin)
 if strcmp(T,'usda')==0
    error('only usda defined so far...')
 end
 if isempty(sand)|isnan(sand)|min(sand)<0
     error('wrong input: sand')
 end
 if isempty(clay)|isnan(clay)|min(clay)<0
     error('wrong input: sand')
 end
 if max(sand)>1
     sand=sand/100;
     warning('sand is >1; ...divided by 100 to get fraction instead of %')
 end
 if max(clay)>1
     clay=clay/100;
     warning('clay is >1; ...divided by 100 to get fraction instead of %')
 end
 if max(sand+clay)>1
    error('data inconsisent: (Clay + Sand)>1')
 end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
silt=1-sand-clay;
i=1;
 for i=1:length(clay)
     if (silt(i)+1.5*clay(i))<.15
         SC{i,:}='SAND';
     elseif ((silt(i)+1.5*clay(i))>=.15)&((silt(i)+2*clay(i))<.3)
         SC{i,:}='LOAMY SAND';
     elseif (clay(i)>=0.07) & (clay(i)<=0.2) & (sand(i)>0.52) & ((silt(i)+2*clay(i))>=0.3) 
         SC{i,:}='SANDY LOAM';
     elseif (clay(i)<0.07) & (silt(i) < 0.5) & ((silt(i)+2*clay(i))>=0.3)
         SC{i,:}='SANDY LOAM';
     elseif (clay(i)>=0.07) & (clay(i)<=0.27) & (silt(i)>=0.28) & (silt(i)<0.5) & (sand(i)<=0.52)
         SC{i,:}='LOAM';
     elseif ((silt(i)>=0.5) & clay(i)>=0.12 & clay(i)<0.27) | (silt(i)>=0.5 & silt(i)<0.8 & clay(i)<0.12)
         SC{i,:}='SILT LOAM';
     elseif silt(i)>=0.8 & clay(i)<0.12
         SC{i,:}='SILT';
     elseif clay(i)>=0.2 & clay(i)<0.35 & silt(i)<0.28 & sand(i)>0.45
         SC{i,:}='SANDY CLAY LOAM';
     elseif clay(i)>=0.27 & clay(i) <0.4 & sand(i)>0.2 & sand(i)<=0.45
         SC{i,:}='CLAY LOAM';
     elseif clay(i)>=0.27 & clay(i)<0.4 & sand(i)<=0.2
         SC{i,:}='SILTY CLAY LOAM';
     elseif clay(i)>=0.35 & sand(i)>=0.45
         SC{i,:}='SANDY CLAY';
     elseif clay(i)>=0.4 & silt(i)>=0.4
         SC{i,:}='SILTY CLAY';
     elseif clay(i)>= 0.4 & sand(i)<=0.45 & silt(i)<0.4
         SC{i,:}='CLAY';
     else
         warning('no soil class found')
     end
 end
 %------------------------------------------------------------------------
 if varargin{1}==1
     close
%      f1=figure;
     f1 = figure(...
    'Color',[1,1,1],...
    'Units','centimeters',... % centimeters, pixels, points, characters
    'Position',[0,0,8,8],...%[left bottom width height]
    'PaperSize',[8,8]);
     axes1 = axes('parent',f1);
     %set(gcf,'Color','w','position',[2800 300 600 600])
     plot([0 100],[0 0],'k');hold on
     set(gca,'XDir','reverse')
     plot([0 50],[0 100],'k')
     plot([50 100],[100 0],'k')
     xlabel('Sand [%]')
     sizelabels = 6;
     colorlabels =[0.1,0.1,0.1];
%      text(95,60,'Clay [%]','FontSize',10)
%      text(25,60,'Silt [%]','FontSize',10)
%      text(55,-5,'Sand [%]','FontSize',10)
%      p1=fill([100 85 95 100],[0 0 10 0],[1 1 1])
%      plot([100 85 95 100],[0 0 10 0],'k','LineWidth',2)
     plot([85 95],[0 10],'k','LineWidth',2)
     %p1=polyshape([100 85 95 100],[0 0 10 0],[1 1 1])
     text(97,3, 'SAND','color',colorlabels,'FontSize',sizelabels)
     
%      p2= fill([85 70 92.4 95 85],[0 0 15 10 0],[1 1 1])
%      plot([85 70 92.4 95 85],[0 0 15 10 0],'k','LineWidth',2)
     plot([70 92.4],[0 15],'k','LineWidth',2)
     text(90,4.5,{'LOAMY';'      SAND'},'color',colorlabels,'FontSize',sizelabels)
     
%      p3=fill([70 50 46.7 55 62 90 92.4 70],[0 0 7 7 20 20 15 0],[1 1 1])
%      plot([70 50 46.7 55 62 90 92.4 70],[0 0 7 7 20 20 15 0],'k','LineWidth',2)
     plot([50 46.7 55 62 90],[0 7 7 20 20],'k','LineWidth',2)
     text(77,10,'SANDY LOAM','color',colorlabels,'FontSize',sizelabels)  
     
%      p4=fill([55 46.7 38 58.1 62 55],[7 7 27 27 20 7],[1 1 1])
%      plot([55 46.7 38 58.1 62 55],[7 7 27 27 20 7],'k')
     plot([46.7 38 58.1 62],[7 27 27 20],'k','LineWidth',2)
     text(55,18,'LOAM','color',colorlabels,'FontSize',sizelabels)
         
%      p5=fill([50 20 13.5 6 13.2 38 50],[0 0 12 12 27 27 0],[1 1 1])
%      plot([50 20 13.5 6 13.2 38 50],[0 0 12 12 27 27 0],'k','LineWidth',2)
     plot([20 13.5 6],[0 12 12],'k','LineWidth',2)
     plot([13.2 38],[27 27],'k','LineWidth',2)
     text(35, 11, 'SILT LOAM','color',colorlabels,'FontSize',sizelabels)
%      p6=fill([20 0 6 13.5 20],[0 0 12 12 0],[1 1 1])
%      plot([20 0 6 13.5 20],[0 0 12 12 0],'k','LineWidth',2)

     text(12,5,'SILT','color',colorlabels,'FontSize',sizelabels)
     
%      p7=fill([90 62 58.1 62.5 82.5 90],[20 20 27 35 35 20],[1 1 1])
%      plot([90 62 58.1 62.5 82.5 90],[20 20 27 35 35 20],'k','LineWidth',2)
     plot([58.1 62.5 82.5],[27 35 35],'k','LineWidth',2)
     text(85,27,'SANDY CLAY LOAM','color',colorlabels,'FontSize',sizelabels)
%      p8=fill([58.1 33 40 65 58],[27 27 40 40 27],[1 1 1])
%      plot([58.1 33 40 65 58],[27 27 40 40 27],'k','LineWidth',2)
     plot([33 40 65 58],[27 40 40 27],'k','LineWidth',2)
     text(55,34,'CLAY LOAM','Color','w','color',colorlabels,'FontSize',sizelabels)
%      p9=fill([33 13.5 20 40 33],[27 27 40 40 27],[1 1 1])
%      plot([33 13.5 20 40 33],[27 27 40 40 27],'k','LineWidth',2)
     plot([20 40],[40 40],'k','LineWidth',2)
     text(33,33,{'     SILTY';'CLAY LOAM'},'color',colorlabels,'FontSize',sizelabels)
     
%      p10=fill([82.5 62.5 72.5 82.5],[35 35 55 35],[1 1 1])
%      plot([82.5 62.5 72.5 82.5],[35 35 55 35],'k','LineWidth',2)
     plot([62.5 72.5],[35 55],'k','LineWidth',2)
     text(77,40,{'SANDY'; ' CLAY'},'Color','w','color',colorlabels,'FontSize',sizelabels)
     
%      p11=fill([65 40 30 50 72.5 65],[40 40 60 100 55 40],[1 1 1])
     plot([40 30],[40 60],'k','LineWidth',2)
     text(55,60,'CLAY','color',colorlabels,'FontSize',sizelabels)
     
%      p12=fill([40 20 30 40],[40 40 60 40],[1 1 1])
%      plot([40 20 30 40],[40 40 60 40],'k')
     text(33,47,{'SILTY' ;'CLAY'},'color',colorlabels,'FontSize',sizelabels)
     %axis square
     daspect([0.866 1 1])
     axis off
     ax = gca;
     ax.InnerPosition=[0.02 0.02 0.96 0.96];
%      -0.167741935483871,-0.134969325153374,1.290322580645161,1.226993865030675
     % Create text
    text('Parent',axes1,'Rotation',60,'String','CLAY (%)',...
    'Position',[83 42 0],...
    'Visible','on','FontSize',8);

% Create text
    text('Parent',axes1,'Rotation',-60,'String','SILT (%)',...
    'Position',[22 53 0],...
    'Visible','on','FontSize',8);

% Create text
    text('Parent',axes1,'String','SAND (%)',...
    'Position',[56 -3 0],...
    'Visible','on','FontSize',8);
     
     i=1;
     for i=1:length(sand)
         o=100-(1-clay(i))*100; %offset
         plot(sand(i)*100+o/2,clay(i)*100,'o','MarkerFaceColor','none','MarkerEdgeColor','k','MarkerSize',3)
     end
 end
 
 %% get SCINT
 SCINT = nan(size(SC,1),1);
 SCINT(strcmp(SC,'SAND')) = 1;
 SCINT(strcmp(SC,'LOAMY SAND')) = 2;
 SCINT(strcmp(SC,'SANDY LOAM')) = 3;
 SCINT(strcmp(SC,'LOAM')) = 4;        
 SCINT(strcmp(SC,'SILT LOAM')) = 5;   
 SCINT(strcmp(SC,'SILT')) = 6;   
 SCINT(strcmp(SC,'SANDY CLAY LOAM')) = 7;   
 SCINT(strcmp(SC,'CLAY LOAM')) = 8;   
 SCINT(strcmp(SC,'SILTY CLAY LOAM')) = 9;   
 SCINT(strcmp(SC,'SANDY CLAY')) = 10;   
 SCINT(strcmp(SC,'SILTY CLAY')) = 11;   
 SCINT(strcmp(SC,'CLAY')) = 12;   
 
 %% get SC2:
 SC2 = SC;
 SC2 = strrep(SC2,'SAND','S');
 SC2 = strrep(SC2,'LOAMY SAND','LS');
 SC2 = strrep(SC2,'SANDY LOAM','SL');
 SC2 = strrep(SC2,'LOAM','L');
 SC2 = strrep(SC2,'SILT LOAM','SIL');
 SC2 = strrep(SC2,'SILT','SI');
 SC2 = strrep(SC2,'SANDY CLAY LOAM','SCL');
 SC2 = strrep(SC2,'CLAY LOAM','CL');
 SC2 = strrep(SC2,'SILTY CLAY LOAM','SICL');
 SC2 = strrep(SC2,'SANDY CLAY','SC');
 SC2 = strrep(SC2,'SILTY CLAY','SIC');
 SC2 = strrep(SC2,'CLAY','C');
 SC2 = strrep(SC2,' ','');
 SC2 = strrep(SC2,'Y','');
 %test: [SC,~,~,SC2] = soil_classification([0:0.1:1,0:0.1:0.5],[1:-0.1:0,0:0.1:0.5],'usda',0)
end %of function

