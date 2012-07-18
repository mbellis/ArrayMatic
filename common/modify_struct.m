function NewStruct=modify_struct(Struct,FirstIndexFlag)
if FirstIndexFlag
    %transfort a like Struct.field(k) into a sructure like Struct(k).field
else
    %transfort a like Struct(k).field into a sructure like Struct(.field(k)
    FieldNames=fields(Struct(1));
    NewStruct=[];
    IsNumber=zeros(size(FieldNames));
    for FieldL=1:length(FieldNames)
        eval(sprintf('Value=Struct(1).%s;',FieldNames{FieldL}));
        if isnumeric(Value)
            IsNumber(FieldL)=1;
        end
    end
    for ItemL=1:length(Struct)
        for FieldL=1:length(FieldNames)
            if IsNumber(FieldL)
                try
                eval(sprintf('NewStruct.%s(%u)=Struct(%u).%s;',FieldNames{FieldL},ItemL,ItemL,FieldNames{FieldL}));
                catch
                    'stop'
                end
            else
                try
                eval(sprintf('NewStruct.%s{%u}=Struct(%u).%s;',FieldNames{FieldL},ItemL,ItemL,FieldNames{FieldL}));
                catch
                    'stop'
                end
            end
        end
    end
end
