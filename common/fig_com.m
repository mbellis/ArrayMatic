% FIG_COM
% functions for image and plot set up

function varargout=fig_com(Action,varargin)
global P
switch Action
    case 'legend'
        h=varargin{1};
        if nargin==3
            Location=varargin{2};
        else
            Location='Best';
        end
        figure(h)
        if P.flag.testAlgo
            Algos=unique(P.point.algo);
            AlgoNb=length(Algos);
            Colors=colors(colormap,AlgoNb);
            Legend={P.point.algo{1}};
            plot(0,0,'color',Colors(1,:))
            AlgoRank=1;
            for PointL=1:P.point.nb
                if isempty(strmatch(P.point.algo{PointL},Legend,'exact'))
                    Legend{end+1}=P.point.algo{PointL};
                    AlgoRank=AlgoRank+1;
                    plot(0,0,'color',Colors(AlgoRank,:))
                end
            end
            legend(Legend,'location',Location)
        else
            BiolNb=P.biol.nb;
            Colors=colors(colormap,BiolNb);
            Legend={strrep(P.biol.name{P.point.biolRank(1)},'_',' ')};
            plot(0,0,'color',Colors(1,:))
            BiolRank=1;
            for PointL=1:P.point.nb
                if isempty(strmatch(strrep(P.biol.name{P.point.biolRank(PointL)},'_',' '),Legend,'exact'))
                    Legend{end+1}=strrep(P.biol.name{P.point.biolRank(PointL)},'_',' ');
                    BiolRank=BiolRank+1;
                    plot(0,0,'color',Colors(BiolRank,:))
                end
            end
            legend(Legend,'location',Location)
        end
        varargout{1}=Legend;
        varargout{2}=Colors;
end