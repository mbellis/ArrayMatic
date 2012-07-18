% COLORS - construct a ItemNbx3 color map by sampling an existing colormap

% c) Michel Bellis
% arraymatic@gmail.com

% NO INPUT
%ColorMap: a ColorNbx3 matrix
%ItemNb: the size of the sampled color map

% OUTPUT
%Colors: the sampled color map

% VERSIONS
%
% V01 - 2010 05 14

function [Colors]=colors(ColorMap,ItemNb)
ColorNb=size(ColorMap,1);
while ItemNb>ColorNb
    if ItemNb-ColorMap>ColorMap
        ColorMap=[ColorMap;ColorMap];
    else
        ColorMap=[ColorMap;ColorMap(1:ItemNb-ColorMap,:)];
    end
    ColorNb=size(ColorMap,1);
end    
ColorStep=ColorNb/max(1,(ItemNb-1));
ColorPos=round([1:ColorStep:ColorNb]);
if length(ColorPos)<ItemNb
    ColorPos(end+1)=ColorNb;
end
Colors=ColorMap(ColorPos,:);

%ColorStep=ColorNb/ItemNb;
% Colors=zeros(ItemNb,3);
% for ItemL=1:ItemNb
%     if mod(round((ItemL-1)*ColorStep)+1,ColorNb)
%         Colors(ItemL,:)=ColorMap(mod(round((ItemL-1)*ColorStep)+1,ColorNb),:);
%     else
%         Colors(ItemL,:)=ColorMap(end,:);
%     end
% end
