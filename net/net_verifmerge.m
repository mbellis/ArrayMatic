%% m95 et m6 et m92
figure
set(gcf,'color',[1,1,1])

RowPsRank=1:200;
ColPsRank=1:200;
C=load_cvm(95,1,RowPsRank,ColPsRank,1,1,1);
subplot(2,4,1)
image(C)
title('m95n1')
C=load_cvm(95,2,RowPsRank,ColPsRank,1,1,1);
subplot(2,4,2)
image(C)
title('m95n2')
C=load_cvm(95,3,RowPsRank,ColPsRank,1,1,1);
subplot(2,4,3)
image(C)
title('m95n3')
C=load_cvm(95,4,RowPsRank,ColPsRank,1,1,1);
subplot(2,4,4)
image(C)
title('m95n4')

ChipRanks=[92,92,92,6];
NetRanks=[1,2,3,63];
ColPos=[1,1,1,2];
MergeFile={'m92n4_m6n63_combinedps_corr60','m92n4_m6n63_combinedps_corr60','m92n4_m6n63_combinedps_corr60','m92n4_m6n63_combinedps_corr60'};
for ChipL=1:length(ChipRanks)
    cd(K.dir.chip)
    eval(sprintf('load %s',MergeFile{ChipL}))
    C=load_cvm(ChipRanks(ChipL),NetRanks(ChipL),PsRanks(RowPsRank,ColPos(ChipL)),PsRanks(ColPsRank,ColPos(ChipL)),1,1,1);
    subplot(2,4,ChipL+4)
    image(C)
    title(sprintf('m%un%u',ChipRanks(ChipL),NetRanks(ChipL)))
end
%% m94 et m6 et m93
figure
set(gcf,'color',[1,1,1])

RowPsRank=1:200;
ColPsRank=1:200;
C=load_cvm(94,1,RowPsRank,ColPsRank,1,1,1);
subplot(2,3,1)
image(C)
title('m94n1')
C=load_cvm(94,2,RowPsRank,ColPsRank,1,1,1);
subplot(2,3,2)
image(C)
title('m94n2')
C=load_cvm(94,3,RowPsRank,ColPsRank,1,1,1);
subplot(2,3,3)
image(C)
title('m94n3')

ChipRanks=[93,93,6];
NetRanks=[1,2,63];
ColPos=[1,1,2]
MergeFile={'m93n3_m6n63_combinedps_corr60','m93n3_m6n63_combinedps_corr60','m6n63_m93n3_combinedps_corr60'};
for ChipL=1:length(ChipRanks)
    cd(K.dir.chip)
    eval(sprintf('load %s',MergeFile{ChipL}))
    C=load_cvm(ChipRanks(ChipL),NetRanks(ChipL),PsRanks(RowPsRank,1),PsRanks(ColPsRank,1),1,1,1);
    subplot(2,3,ChipL+3)
    image(C)
    title(sprintf('m%un%u',ChipRanks(ChipL),NetRanks(ChipL)))
end




%% m93 et m91 et m92

figure
set(gcf,'color',[1,1,1])

RowPsRank=1:200;
ColPsRank=1:200;
C=load_cvm(93,1,RowPsRank,ColPsRank,1,1,1);
subplot(2,2,1)
image(C)
title('m93n1')
C=load_cvm(93,2,RowPsRank,ColPsRank,1,1,1);
subplot(2,2,2)
image(C)
title('m93n2')

ChipRanks=[91,92];
NetRanks=[3,4];
ColPos=[1,2];
MergeFile={'m91n3_m92n4_combinedps_corr60','m91n3_m92n4_combinedps_corr60'};
for ChipL=1:length(ChipRanks)
    cd(K.dir.chip)
    eval(sprintf('load %s',MergeFile{ChipL}))
    C=load_cvm(ChipRanks(ChipL),NetRanks(ChipL),PsRanks(RowPsRank,ColPos(ChipL)),PsRanks(ColPsRank,ColPos(ChipL)),1,1,1);
    subplot(2,2,ChipL+2)
    image(C)
    title(sprintf('m%un%u',ChipRanks(ChipL),NetRanks(ChipL)))
end


%% m92 et m5 m8 et m27

figure
set(gcf,'color',[1,1,1])

RowPsRank=1:200;
ColPsRank=1:200;
C=load_cvm(92,1,RowPsRank,ColPsRank,1,1,1);
subplot(2,3,1)
image(C)
title('m92n1')
C=load_cvm(92,2,RowPsRank,ColPsRank,1,1,1);
subplot(2,3,2)
image(C)
title('m92n2')
C=load_cvm(92,3,RowPsRank,ColPsRank,1,1,1);
subplot(2,3,3)
image(C)
title('m92n3')

ChipRanks=[5,8,27];
NetRanks=[123,228,164];
ColPos=[1,1,1];
MergeFile={'m5_combinedps_corr60','m8_combinedps_corr60','m27_combinedps_corr60'};
for ChipL=1:length(ChipRanks)
    cd(K.dir.chip)
    eval(sprintf('load %s',MergeFile{ChipL}))
    C=load_cvm(ChipRanks(ChipL),NetRanks(ChipL),PsRanks(RowPsRank,ColPos(ChipL)),PsRanks(ColPsRank,ColPos(ChipL)),1,1,1);
    subplot(2,3,ChipL+3)
    image(C)
    title(sprintf('m%un%u',ChipRanks(ChipL),NetRanks(ChipL)))
end

%% m91 et m2 et m3

figure
set(gcf,'color',[1,1,1])

RowPsRank=1:200;
ColPsRank=1:200;
C=load_cvm(91,1,RowPsRank,ColPsRank,1,1,1);
subplot(2,3,1)
image(C)
title('m91n1')
C=load_cvm(91,2,RowPsRank,ColPsRank,1,1,1);
subplot(2,3,2)
image(C)
title('m91n2')
C=load_cvm(91,3,RowPsRank,ColPsRank,1,1,1);
subplot(2,3,3)
image(C)
title('m91n3')


ChipRanks=[2,3];
NetRanks=[80,86];
ColPos=[1,1]
MergeFile={'m91n3_m92n4_combinedps_corr60','m91n3_m92n4_combinedps_corr60'};
for ChipL=1:length(ChipRanks)
    cd(K.dir.chip)
    eval(sprintf('load %s',MergeFile{ChipL}))
    C=load_cvm(ChipRanks(ChipL),NetRanks(ChipL),PsRanks(RowPsRank,ColPos(ChipL)),PsRanks(ColPsRank,ColPos(ChipL)),1,1,1);
    subplot(2,3,ChipL+3)
    image(C)
    title(sprintf('m%un%u',ChipRanks(ChipL),NetRanks(ChipL)))
end