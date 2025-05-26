function Positions = DivideScreen(GrR,GrC,ScSz, w_padd,h_padd, SubSc)
% Positions = DivideScreen(GrR,GrC,ScSz, w_padd,h_padd, SubSc)
%     The function will divide the screen to a GrR (rows) by GrC (columns) grid with each
%     figure padded with a horizontal padding (w_padd) and a vertical padding
%     (h_padd). In order to support plotting multiple grids in one figure
%     you can enter a sub area of the screen to be divided (SubSc).
%
%     GrR, GrC - a GrC by GrR grid
%     ScSz     - Screen size in the matlab format (i.e. [left bottom width height])
%     w_padd   - figure padding (left+right padding)
%     h_padd   - figure padding (bottom+top padding)
%     SubSc    - Sub area of the screen to be divided in the matlab format (i.e. [left bottom width height])
%     (all units are in pixels)
%
%
%     Positions - output cell aray of size (1,GrR*GrC) and indexing
%                 compatible with matlab's subplot command
%
%
% Examples:
%
%   figure; set(gcf,'Position',get(0,'ScreenSize'));
%   Positions = DivideScreen(4,3,get(0,'ScreenSize'),30,30);
%   for i = 1:12
%         subplot('Position',Positions{i});
%   end
%
%
%  figure; set(gcf,'Position',get(0,'ScreenSize'));
%  Positions = DivideScreen(3,2,get(0,'ScreenSize'),10,10,[100 100 300 300]);
%  for i = 1:6
%        subplot('Position',Positions{i});
%  end
%
%
%
if nargin<6
    StartX = ScSz(1);
    StartY = ScSz(2);
    lX   = ScSz(3);
    lY   = ScSz(4);
else
    StartX = SubSc(1);
    StartY = SubSc(2);
    lX   = SubSc(3);
    lY   = SubSc(4);
end
if nargin<4
    w_padd = lX/GrC/10;
    h_padd = lY/GrR/10;
end
w_padd = w_padd/ScSz(3);
h_padd = h_padd/ScSz(3);
lX   = lX/(ScSz(3))-w_padd/2;
lY   = lY/(ScSz(4))-h_padd/2;
StartX = StartX/(ScSz(3))+w_padd/2;
StartY = StartY/(ScSz(4))+h_padd/2;
dX = (lX)/GrC;
dY = (lY)/GrR;

fig_w  = dX - w_padd;
fig_h  = dY - h_padd;
shiftX = w_padd+fig_w;
shiftY = h_padd+fig_h;

Positions = cell(1,GrC*GrR);
cnt = GrR;
for i =1:GrR
    for j = 1:GrC
        Positions{(cnt-1)*(GrC)+j} = [StartX+(j-1)*shiftX StartY+(i-1)*shiftY fig_w fig_h];
    end
    cnt = cnt -1;
end

    
    