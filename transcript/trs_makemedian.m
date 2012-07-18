% TRS_MAKEMEDIAN - construct a new point with median rank values

% c) Michel Bellis
% arraymatic@gmail.com

% NO INPUT

% NO OUTPUT

function trs_makemedian()
global P DataRanks
% initialise VarAll
VarAll=[];

%recover ranks from all points grouped by biological conditions
for BiolL=1:P.biol.nb
    Ranks=[];
    CurrGrp=P.biol.pointIndex{BiolL};
    PointNb=length(CurrGrp);
    if PointNb>0
        for PointL=1:PointNb
            CurrRank=CurrGrp(PointL);
            if P.flag.loadData
                CurrRanks=load_data('DataRanks.float32le',P.dir.data,P.chip.currProbeSetNb,P.point.nb,'single','ieee-le',1:P.chip.currProbeSetNb,CurrRank);
                Ranks=[Ranks;CurrRanks'];                
            else
                Ranks=[Ranks;DataRanks(:,CurrRank)'];
            end
        end
    end

    % increment VarAll with the median rank values of the current
    % biological condition
    if ~isempty(Ranks)
        % not measured values have been replaced by -2 (absent), -1 (not measured), 0 (under a threshold) or
        % 100 (saturated values)
        MissIndex=find(Ranks==100|Ranks==0|Ranks==-1|Ranks==-2);
        if ~isempty(MissIndex)
            Ranks(MissIndex)=NaN;
        end
        if size(Ranks,1)==1
            VarAll=[VarAll;Ranks];
        else
            VarAll=[VarAll;nanmedian(Ranks)];
        end
    end

end


if ~isempty(VarAll)
    %calculate median point
    if size(VarAll,1)==1
        MedianRanks=VarAll;
    else
        MedianRanks=nanmedian(VarAll);
    end

    %convert NaN into -1
    NaNIndex=find(isnan(MedianRanks));
    if ~isempty(NaNIndex)
        MedianRanks(NaNIndex)=-1;
    end
    
    %column vector
    MedianRanks=MedianRanks';

    %save median point
    cd(P.dir.data)
    save MedianRanks MedianRanks    
else
    errordlg('Median Point is empty. Process canceled')
    error('process canceled')
end
