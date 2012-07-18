function modify_s(ResRank)
global S
MemS=S;
S=cell(size(MemS));
for PosL=1:length(MemS)
    S{PosL}=cell(size(MemS{PosL}));
end

%fill S
for PosL1=1:length(MemS)
    if ~isempty(MemS{PosL1})
        for PosL2=1:length(MemS{PosL1})
            if ~isempty(MemS{PosL1}{PosL2})
                CurrMemS=MemS{PosL1}{PosL2};
                CurrS.position=[];
                CurrS.clearFlag=[];
                CurrS.testCurve={};
                CurrS.testSurf=[];
                CurrS.mean={};
                CurrS.std={};
                CurrS.cdf={};
                CurrS.minZVar=[];
                CurrS.maxZVar=[];
                CurrS.zval={};
                Pos=0;
                for LineL=1:size(CurrMemS.testSurf,1)
                    Columns=find(CurrMemS.testSurf(LineL,:));
                    if ~isempty(Columns)
                        for ColumnL=1:length(Columns)
                            Pos=Pos+1;
                            CurrS.position(Pos,:)=[LineL,Columns(ColumnL),Pos];
                            CurrS.clearFlag(Pos,1)=CurrMemS.clearFlag(LineL,Columns(ColumnL));
                            CurrS.testCurve{Pos,1}=CurrMemS.testCurve{LineL,Columns(ColumnL)};
                            CurrS.testSurf(Pos,1)=CurrMemS.testSurf(LineL,Columns(ColumnL));
                            CurrS.mean{Pos,1}=CurrMemS.mean{LineL,Columns(ColumnL)};
                            CurrS.std{Pos,1}=CurrMemS.std{LineL,Columns(ColumnL)};
                            CurrS.cdf{Pos,1}=CurrMemS.cdf{LineL,Columns(ColumnL)};
                            CurrS.minZVar(Pos,1)=CurrMemS.minZVar(LineL,Columns(ColumnL));
                            CurrS.maxZVar(Pos,1)=CurrMemS.maxZVar(LineL,Columns(ColumnL));
                            CurrS.zval{Pos,1}=CurrMemS.zval{LineL,Columns(ColumnL)};
                        end
                    end
                end
                S{PosL1}{PosL2}=CurrS;
            end
        end
    end
end

cd(P.dir.data)
eval(sprintf('save CalibSet_%02u S',ResRank))
                    


