function Structure=del_structitems(Structure,Pos)
FieldNames=fieldnames(Structure);
Pos=sort(Pos);
for PosL=length(Pos):-1:1
    for FieldL=1:length(FieldNames)
        eval(sprintf('Structure.%s(%u)=[];',FieldNames{FieldL},Pos(PosL)));
    end
end
