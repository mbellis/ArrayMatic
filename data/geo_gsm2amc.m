% recover information used to fill tabular file used to describe a particular GSE
Gpl='GPL1355';
Gse='GSE29488';
cd('E:\sosma\collab\2012_scamps');
n=textread(sprintf('%s_%s.txt',Gse,Gpl),'%s');
itemnb=length(n);
switch Gse
    case 'GSE29488'
        kin=cell(itemnb,1);
        for i=1:itemnb
            tval=regexp(n{i},'\d+(?=_day)','match');
            if ~isempty(tval)
                kin{i}=[tval{1},'d'];
            else
                kin{i}='-';
            end
        end
        rep=cell(itemnb,1);
        for i=1:itemnb
            tval=regexp(n{i},'(?<=rep)\d+','match');
            if ~isempty(tval)
                rep{i}=[tval{1}];
            else
                rep{i}='-';
            end
        end
        level=cell(itemnb,1);
        for i=1:itemnb
            if ~isempty(findstr('(R)',n{i}))
                level{i}='rostral';
            elseif ~isempty(findstr('(M)',n{i}))
                level{i}='idem';
            elseif ~isempty(findstr('(C)',n{i}))
                level{i}='caudal';
            else
                level{i}='-';
            end
        end
        treat=cell(itemnb,1);
        for i=1:itemnb
            if ~isempty(findstr('(K)',n{i}))
                treat{i}='tube';
            elseif ~isempty(findstr('(S)',n{i}))
                treat{i}='NT3';
            elseif ~isempty(findstr('D)',n{i}))
                treat{i}='cont';
            else
                treat{i}='-';
            end
        end

        biol=cell(itemnb,1);
        for i=1:itemnb
            biol{i}=sprintf('SC_%s_%s_%s',treat{i},level{i},kin{i});
        end
        Biol=unique(biol);

        point=cell(itemnb,1);
        for i=1:itemnb
            point{i}=sprintf('%s_%s',biol{i},rep{i});
        end

    case 'GSE34000'
        kin=cell(itemnb,1);
        for i=1:itemnb
            tval=regexp(n{i},'\d_week','match');
            if ~isempty(tval)
                kin{i}=[tval{1}(1),'w'];
            end
        end
        inj=cell(itemnb,1);
        for i=1:itemnb
            if ~isempty(findstr('Normal',n{i}))
                inj{i}='cont';
            else
                inj{i}='STZ';
            end
        end
        treat=cell(itemnb,1);
        for i=1:itemnb
            if ~isempty(findstr('FK',n{i}))
                treat{i}='FK';
            else
                treat{i}='vehicle';
            end
        end
        biol=cell(itemnb,1);
        for i=1:itemnb
            biol{i}=sprintf('%s_%s_%s',inj{i},treat{i},kin{i});
        end
        Biol=unique(biol);

        rep=zeros(itemnb,1);
        for i=1:length(Biol)
            pos=strmatch(Biol{i},biol,'exact');
            for j=1:length(pos)
                rep(pos(j))=j;
            end
        end

        point=cell(itemnb,1);
        for i=1:itemnb
            point{i}=sprintf('%s_%u',biol{i},rep(i));
        end


    case 'GSE464'
        inj=cell(itemnb,1);
        for i=1:itemnb
            if ~isempty(findstr('_Ctr',n{i}))|~isempty(findstr('-con',n{i}))
                inj{i}='cont';
            elseif ~isempty(regexp(n{i},'sev[_-]','match'))
                inj{i}='sev';
            elseif ~isempty(regexp(n{i},'mod[_-]','match'))
                inj{i}='mod';
            else
                inj{i}='mild';
            end
        end
        loc=cell(itemnb,1);
        for i=1:itemnb
            if ~isempty(regexp(n{i},'InjA|CtrA','match'))
                loc{i}='T8';
            elseif ~isempty(regexp(n{i},'InjB|CtrB','match'))
                loc{i}='T10';
            else
                loc{i}='T9';
            end
        end
        kin=cell(itemnb,1);
        for i=1:itemnb
            tval=regexp(n{i},'\d+d','match');
            if ~isempty(tval)
                kin{i}=tval{1};
            else
                tval=regexp(n{i},'\d+h','match');
                if ~isempty(tval)
                    kin{i}=tval{1};
                else
                    tval=regexp(n{i},'\d+m','match');
                    if ~isempty(tval)
                        kin{i}=tval{1};
                    else
                        kin{i}='t0';
                    end
                end
            end
        end

        biol=cell(itemnb,1);
        for i=1:itemnb
            biol{i}=sprintf('%s_%s_%s',inj{i},loc{i},kin{i});
        end
        Biol=unique(biol);

        rep=zeros(itemnb,1);
        for i=1:length(Biol)
            pos=strmatch(Biol{i},biol,'exact');
            for j=1:length(pos)
                rep(pos(j))=j;
            end
        end

        point=cell(itemnb,1);
        for i=1:itemnb
            point{i}=sprintf('%s_%u',biol{i},rep(i));
        end

        biolr=zeros(itemnb,1);
        sinj={'cont','mild','mod','sev'};
        sloc={'T8','T9','T10'};
        skin={'t0','30m','4h','12h','24h','48h','72h','4d','7d','14d','21d','28d'};
        biolrank=0;
        for i=1:length(sinj)
            pos1=strmatch(sinj{i},inj,'exact');
            if~isempty(pos1)
                for j=1:length(sloc)
                    pos2=strmatch(sloc{j},loc,'exact');
                    if~isempty(intersect(pos1,pos2))
                        for k=1:length(skin)
                            pos3=strmatch(skin{k},kin,'exact');
                            pos=intersect(pos1,intersect(pos2,pos3));
                            if~isempty(pos)
                                biolrank=biolrank+1;
                                biolr(pos)=biolrank;
                            end
                        end
                    end
                end
            end
        end
end