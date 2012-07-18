function set_figsize(Type)
% global P
% h=P.handle
set(gcf,'units','pixels')
switch Type
    case '1024px'        
        set(gcf,'position',[1,1,1020,764])
    case '1280px'
        set(gcf,'position',[1,1,1280,956])
    case '1910px'
        set(gcf,'position',[1,1,1910,1064])
    case '780px'
        set(gcf,'position',[1,1,780,780])
    case '960px'
        set(gcf,'position',[1,1,960,960])

end