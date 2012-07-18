% SELECT_NET allows to select one or several existing networks

% INPUT
% SelectType: either 'unique' or 'multiple'
% varargin: 
% 1 ModelRank: rank of the chip set model

% OUPUT
% ModelRank: rank of the chip set model
% NetRank: rank(s) of the selected net(s) corresponding to the net(s) name
%          (sprintf('n%05u',NetRank)
% NetPos: the positions(s) of the selected net(s) in K.net{ModelRank}

function [ModelRank,NetRank,NetPos]=select_net(SelType,varargin)
global K

if ~isequal(SelType,'unique')&&~isequal(SelType,'multiple')
    h=errodlg(sprintf('SelType must be unique or multiple. %s not allowed',SelType));
    waitfor(h)
    error('process canceled')
end
if nargin==1
[ModelRank,Ok]=listdlg('liststring',K.chipSet.name,'selectionmode','single','listsize',[400 300],'promptstring','Select chip model','name','COM_RAY');
else
    ModelRank=varargin{1};
    Ok=1;
end
if Ok==1
    NetMade=find(K.net{ModelRank}.netMade==1);
    NetRanks=K.net{ModelRank}.rank(NetMade);
    if ~isempty(NetRanks)
        NetList=K.net{ModelRank}.name(NetMade);
        [SelNetPos,Ok]=listdlg('liststring',NetList,'selectionmode',SelType,'listsize',[400 300],'promptstring','Select net','name','COM_RAY');
        if Ok==1
            NetNb=length(SelNetPos);
            NetRank=zeros(NetNb,1);
            NetPos=zeros(NetNb,1);
            for NetL=1:NetNb
                NetRank(NetL)=NetRanks(SelNetPos(NetL));
                NetPos(NetL)=find(K.net{ModelRank}.rank==NetRank(NetL));
            end
        end
    end
else
    ModelRank=0;
    NetRank=0;
    NetPos=0;
end
