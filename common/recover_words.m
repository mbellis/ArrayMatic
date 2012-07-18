%========================%
% FUNCTION RECOVER_WORDS %
%========================%
%
% RECOVER_WORDS recovers all the words in a string. The program is adapted to process GEO 
% records and has been set up after a thorough analysis of its characteristics.
%
%INPUT PARAMETERS
%
% 1    String: string to be parsed
% 2 Separator: separator between words in numeric format (9 = tabulation, 32 = space,
%              59= semicolon)
% 3   MinSize: minimal size of words
% 4   MaxSize: maximal size of words (no MaxSize if =0)
%
%OUTPUT PARAMETERS
%
% 1 WordList: words found in String
% 2   WordNb: occurence of each word
% 3 WordSize: length of words


%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%
%                          c) Michel Bellis                                                %
%                          michel.bellis@crbm.cnrs.fr                                      %
%            Affiliation:  CNRS (Centre National de la Recherche Scientifique - France)    %
%  Bioinformatic Project:  ARRAYMATIC => http://code.google.com/p/arraymatic               %
%        Code Repository:  GITHUB => http://github.com/mbellis                             %
%          Personal Page:  http://bns.crbm.cnrs.fr                                         %
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%  THIS CODE IS DISTRIBUTED UNDER THE CeCILL LICENSE, WHICH IS COMPATIBLE WITH       %
%  THE GNU GENERAL PUBLIC LICENCE AND IN ACCORDANCE WITH THE EUROPEAN LEGISLATION.   %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

function [WordList,WordNb,WordSize]=recover_words(String,Separator,MinSize,MaxSize)

%find literal form of the separator
switch Separator
    case 8
        Separator='\b';
    case 9
        Separator='\t';
    case 10
        Separator='\n';
    case 10
        Separator='\v';        
    case 12
        Separator='\f';
    case 13
        Separator='\r';
    otherwise
        Separator=regexptranslate('escape',char(Separator));
end

% activate to stop on string where a problem has been detected
% if ~isempty(findstr(String,'1-(2,4-dichlorobenzyl)-1H-indazole-3-carbohydrazide'))|~isempty(findstr(String,'1,2-dichlorobenzene'))|~isempty(findstr(String,'3H-1,2-dithiole-3-thione'))
%     'stop'
% end

WordList={};
WordNb=[];
WordSize=[];

%replace series of non printable characters, tabulation, space...(char(1) to char(32))
%by separator 
String=regexprep(String,'[\x1-\x{20}]+',Separator);
String=regexprep(String,sprintf('%s+',Separator),Separator);
% " replaced by separator
String=regexprep(String,'"',sprintf('%s',Separator));
% {} (replacement to keep ER-alpha ou IFN-gamma)
String=regexprep(String,'[{}]+','');
% ! replaced by separator
String=regexprep(String,'!',sprintf('%s',Separator));
% $ replaced by ? (exist TM$-ERT cell line)
String=regexprep(String,'\$+','?');
% £ replaced by separator (exist FDR£0.01)
String=regexprep(String,'£+',sprintf('%s',Separator));
% § replaced by nothing 
String=regexprep(String,'§+','');
%final period deleted
String=regexprep(String,sprintf('\.%s*$',Separator),'');
% eliminate separator(s) after 3' and 5'
String=regexprep(String,sprintf('(?<=[3,5]'')%s+',Separator),'');
% replace char(8722) (-)
String=regexprep(String,'\x2212','-');

%%process Keywords
% [ExpStart,EndExp]=regexp(String,sprintf('[Kk]eywords?%s*[=:]%s*',Separator,Separator),'start','end');
% if ~isempty(EndExp)
%     CurrString=String(EndExp(1):end);
%     CurrString=regexprep(CurrString,sprintf(';*%s*(([Kk]eywords?%s*[=:]%s*)|([/;]%s*))',Separator,Separator,Separator,Separator),'*');
%     CurrString=['*',CurrString,'*'];
%     Matches=regexp(CurrString,sprintf('(?<=*)[^*%s]+(?=*)',Separator),'match');
%     for WordL=1:length(Matches)
%         Word=Matches{WordL};
%         Word=regexprep(Word,'^[\\.;:\\]\\[)(,]+','');
%         Word=regexprep(Word,'[\\.;:\\]\\[)(,]+$','');        
%         KeywordList{end+1,1}=Word;
%     end       
%     KeywordList=unique(KeywordList);
%     KeywordNb=ones(length(KeywordList),1);
%     KeywordSize=zeros(length(KeywordList),1);
%     for KeyL=1:length(KeywordList)
%         KeywordSize(KeyL)=length(KeywordList{KeyL});
%     end
%     %delete Keyword part
%     String=String(1:ExpStart(1)-1);
% end

%add a separator at the begining and at the end to process correctly last word
String=regexprep(String,sprintf('^([^%s])',Separator),sprintf('%s$1',Separator));
String=regexprep(String,sprintf('([^%s])$',Separator),sprintf('$1%s',Separator));
%lower eventually first letter
String=regexprep(String,sprintf('^%s([A-Z])',Separator),sprintf('%s${lower($1)}',Separator));

%replace several separators or semi-colons that are tabulation marks (';Separator') and 
% period, comma or colon followed by a Separator by a Separator
% lower eventually first letter
String=regexprep(String,sprintf(',+%s+',Separator),sprintf('%s',Separator));
String=regexprep(String,sprintf('[\\.;:?]+%s+([A-Z]?)',Separator),sprintf('%s${lower($1)}',Separator));

% ~ replaced by nothing if preceded by separator (either followed by a separator or by a
% number)
String=regexprep(String,sprintf('%s+~',Separator),sprintf('%s',Separator));
% # sticked to following numbers 
String=regexprep(String,sprintf('(?<=#)%s+(?=\\d)',Separator),'');
% and deleted in front of a letter
String=regexprep(String,'#+(?=[a-zA-Z])','');
% & replaced by a separator if preceded by separator
String=regexprep(String,sprintf('%s+&',Separator),sprintf('%s',Separator));
% numbers sticked to following °C
String=regexprep(String,sprintf('(?<=\\d)%s+(?=°C)',Separator),'');
% numbers sticked to following %
String=regexprep(String,sprintf('(?<=\\d)%s+(?=%%)',Separator),'');
% % separated from following letter
String=regexprep(String,'%(?=[a-zA-Z])',sprintf('%%%s',Separator));
% 1 x 10^ formated to 1x10^
String=regexprep(String,sprintf('(?<=\\d)%s*[xX]%s*(?=10\\^)',Separator,Separator),'x');
% | replaced by separator if flanked by separators
String=regexprep(String,sprintf('%s+\\|%s+',Separator,Separator),sprintf('%s',Separator));
% process p value => p
String=regexprep(String,'[pP]-?[vV]alue','p');
% process time => t
String=regexprep(String,sprintf('(?<=%s+)[tT]ime(?=%s*=)',Separator,Separator),'t');
%protect Keyword = flank = with separators to find non processed = in the following command
String=regexprep(String,sprintf('[kK]eywords%s*=%s*',Separator,Separator),'keywords=');
String=regexprep(String,'=',sprintf('%s=%s',Separator,Separator));
% =,<> stick n,N,p,P,r,R,r2,R2,t,T and numbers to sign equal,< or > 
String=regexprep(String, sprintf('((%s+([nNpPrRtT]|[rR]2))%s*([=<>]+)%s*(?=\\d)',Separator,Separator,Separator),sprintf('%s${lower($1)}$2',Separator));
% replace non processed = with a separator
String=regexprep(String,sprintf('%s+=%s+',Separator,Separator),sprintf('%s',Separator));
%replace * before or inside numbers by -
String=regexprep(String,'(?<=\d)*(?=\d)','-');
String=regexprep(String,'(?<=-)*(?=\d)','');
%clear all other *
String=regexprep(String,'*','');
%grave accent plus s (error in disease name: e.g. Hodgkin`s lymphoma)
String=regexprep(String,sprintf('`+s%s+',Separator),sprintf('%s',Separator));
%accent plus s (disease name: e.g.  Parkinson's disease)
String=regexprep(String,sprintf('''+s%s+',Separator),sprintf('%s',Separator));
% m2,m3
String=regexprep(String,'²','2');
String=regexprep(String,'([m|cm|mm|µm])\^([2,3])','$1$2');
%replace : (now always in \w:\w expression) by separator
String=regexprep(String,':',sprintf('%s',Separator));

%COMPLEX CHARACTERS
% / 
%recover units (uses star as delimiter since all stars have been cleared)
%three units
String=regexprep(String,sprintf('(?<=\\d%s*)([mµunp]?g|micrograms?)/([kK]g)[/-](min|day|d|week|wk|year|yr|session)',Separator),sprintf('%s$1*$2*$3%s',Separator,Separator));
%two units
String=regexprep(String,sprintf('(?<=\\d%s*)([mµunp]?g|microg\\w*|([m,µ,u](mole|mol)))/(([kK]g)|([dm][lL]))',Separator),sprintf('%s$1*$2%s',Separator,Separator));
String=regexprep(String,sprintf('(?<=\\d%s*)(ml|mL|([m,µ,u](mole|mol))|mOsm)/([kK]g)',Separator),sprintf('%s$1*$2%s',Separator,Separator));
String=regexprep(String,sprintf('(?<=\\d%s*)(I?U)/(ml|g)',Separator),sprintf('%s$1*$2%s',Separator,Separator));
String=regexprep(String,sprintf('(?<=\\d%s*)(m[lLg]|pfu)/(mouse|animal)',Separator),sprintf('%s$1*$2%s',Separator,Separator));
String=regexprep(String,sprintf('(?<=\\d%s*)([mµunp]?g|microg\\w*|[kK][cC]al)/(min|day|d|week|wk|year|yr|session)',Separator),sprintf('%s$1*$2%s',Separator,Separator));
String=regexprep(String,sprintf('(?<=\\d%s*)(min\\w*|hrs|h|day|d|week|wk)/(day|d|week|wk|year|yr|session)',Separator),sprintf('%s$1*$2%s',Separator,Separator));
String=regexprep(String,sprintf('(?<=\\d%s*)([mµunp]?g|microg\\w*|J|W)/((cm|m)[23])',Separator),sprintf('%s$1*$2%s',Separator,Separator));
%clear +/+ ...
%some name must be recovered
String=regexprep(String,'_/','-/');
String=regexprep(String,'/_','/-');
%eliminate genotype marks
String=regexprep(String,'[+-]/[+-]','');
%process common mouse strains (replace / by a star)
String=regexprep(String,'BALB/[cC]','BALB*c');
String=regexprep(String,'C57B[Ll]?/','C57BL*');
String=regexprep(String,'C57/B[Ll]?/','C57BL*');
String=regexprep(String,'DBA/(?=[1,2])','DBA*');
%replace / by separator
String=regexprep(String,sprintf('%s+/%s*',Separator,Separator),sprintf('%s',Separator));
String=regexprep(String,sprintf('%s*/%s+',Separator,Separator),sprintf('%s',Separator));
%recover /
String=regexprep(String,'*','/');


%.
%correct known missing space after .
String=regexprep(String,'\.((Animal)|(Both)|(Cell)|(Comparison)|(Data)|(Early)|(Following)|(Furthermore)|(Gene)|(However)|(In)|(It)|(No)|(On)|(Our)|(Microarray)|(Principal)|(Prior)|(Recentrly)|(Six)|(Somat)|(Some)|(The)|(These)|(This)|(Thus)|(To)|(Using)|(We))',sprintf('%s${lower($1)}',Separator));
%eliminate list of authors name
String=regexprep(String,sprintf('((%sand)?%s([A-Z]\\.)+(%s[A-Z][a-z]+)+[,\\.])+',Separator,Separator,Separator),'');
%eliminate some abreviations
String=regexprep(String,sprintf('%s(M|Ph|Psy)\\.D\\.?',Separator),'');
String=regexprep(String,sprintf('%sPh\\.D\\.',Separator),'');
%eliminate URLs
String=regexprep(String,sprintf('%s(https?)|(ftp)://(\w+)+%s',Separator,Separator),'');

%-
%exists series of '- - - - - - ' or '-----------'
String=regexprep(String,sprintf('(%s?-|%s?-){2,}',Separator,Separator),sprintf('%s',Separator));
String=regexprep(String,sprintf('(%s+-%s+)',Separator,Separator),sprintf('%s',Separator));
%clear minus sign before numbers
String=regexprep(String,sprintf('%s-(\\(?)(?=\\d)',Separator),sprintf('%s$1',Separator));
% and at end of words
String=regexprep(String,sprintf('(?<=[^%s])-+%s',Separator,Separator),sprintf('%s',Separator));

%some expressions to be abbreviated
String=regexprep(String,sprintf('[fF]alse%s*[dD]iscovery%s*[rR]ate',Separator,Separator),sprintf('%sFDR',Separator));

%recover all words
%add a separator at the begining and at the end to process correctly first and last word
String=regexprep(String,sprintf('^([^%s])',Separator),sprintf('%s$1',Separator));
String=regexprep(String,sprintf('([^%s])$',Separator),sprintf('$1%s',Separator));

%SepPos=findstr(String,Separator);
Matches=regexp(String,sprintf('(?<=%s)[^%s]+(?=%s)',Separator,Separator,Separator),'match');
%for WordL=1:length(SepPos)-1
% if ~isempty(findstr(String,'Hdh(Q92'))|~isempty(findstr(String,'Labeled'))|~isempty(findstr(String,'Leukocyte'))|~isempty(findstr(String,'PyMT)F1'))
%     'stop'
% end

for WordL=1:length(Matches)
    %Word=String(SepPos(WordL)+1:SepPos(WordL+1)-1);
    Word=Matches{WordL};
    %     if ~isempty(strfind(Word,'Targeting'))
    %         Word
    %         GPL
    %         GseL
    %         'stop'
    %     end

    %remove punctuation at start and end of words and abnormal parenthesis or brakets
    Word=regexprep(Word,'^[-_'',;:\s\?\.\*\]\)]+','');
    Word=regexprep(Word,'[_'',;:\s\?\.\*\[\(]+$','');
    %remove 's at end of word
    Word=regexprep(Word,'\s+s$','');
    %remove special quotes
    Word=regexprep(Word,sprintf('^[%c%c%c%c]+',char(8216),char(8217),char(8220),char(8221)),'');
    Word=regexprep(Word,sprintf('[%c%c%c%c]+$',char(8216),char(8217),char(8220),char(8221)),'');
    %remove bracketing parenthesis or brakets
    MemLength=length(Word)+1;
    while MemLength~=length(Word) & length(Word)>1
        MemLength=length(Word);
        
        if isequal(Word(1),'(')&isequal(Word(end),')')
            Word=Word(2:end-1);
        end
        if length(Word)>1
            if isequal(Word(1),'[')&isequal(Word(end),']')
                Word=Word(2:end-1);
            end
        end
    end
    MemLength=length(Word)+1;

    while MemLength~=length(Word) & length(Word)>1
        MemLength=length(Word);
        if ~isempty(regexp(Word,'(?<=^\()[^(]+('))|~isempty(regexp(Word,'(?<=^\()[^)]+$'))
            Word=Word(2:end);
        end
        if ~isempty(regexp(Word,'(?<=^\[)[^\[]+['))|~isempty(regexp(Word,'(?<=^\[)[^\[]+$'))
            Word=Word(2:end);
        end

        if ~isempty(regexp(Word,')[^)]+(?=\)$)'))|~isempty(regexp(Word,'[^)]+(?=\)$)'))
            Word=Word(1:end-1);
        end
        if ~isempty(regexp(Word,'][^\]]+(?=\]$)'))|~isempty(regexp(Word,'[^\]]+(?=\]$)'))
            Word=Word(1:end-1);
        end
    end

    %remove punctuation at start and end of words and abnormal parenthesis or brakets
    Word=regexprep(Word,'^[-_'',;:\s\?\.\*\]\)]+','');
    Word=regexprep(Word,'[_'',;:\s\?\.\*\[\(]+$','');
    %remove 's at end of word
    Word=regexprep(Word,'\s+s$','');
    
    % split series of genes
    CurrMatches={};
    if ~isempty(findstr('+',Word))
        if ~isempty(regexp(Word,'[+-]$'))
            %list of genes followed by + or - (CD4+CD8- => CD4+ CD8-)
            Word=regexprep(Word,'([+-])(?=\w)',sprintf('$1%s',Separator));
        else
            %list of genes (CD4+CD8 => CD4 CD8)
            Word=regexprep(Word,'[+-](?=\w)',sprintf('%s',Separator));
        end
        %recover new formed words
        CurrMatches=regexp(Word,sprintf('(?<=%s?)[^%s]+(?=%s?)',Separator,Separator,Separator),'match');
    end
    if length(CurrMatches)>1
        for WordL1=1:length(CurrMatches)
            %Word=String(SepPos(WordL)+1:SepPos(WordL+1)-1);
            CurrWord=CurrMatches{WordL1};
            if ~isempty(CurrWord)
                WordPos=strmatch(CurrWord,WordList,'exact');
                if isempty(WordPos)
                    WordList{end+1,1}=CurrWord;
                    WordNb(end+1,1)=1;
                    WordSize(end+1,1)=length(CurrWord);
                end
            end
        end
    else
        if ~isempty(Word)
            WordPos=strmatch(Word,WordList,'exact');
            if isempty(WordPos)
                WordList{end+1,1}=Word;
                WordNb(end+1,1)=1;
                WordSize(end+1,1)=length(Word);
            end
        end
    end
end
%filter on size
MinPos=find(WordSize<MinSize);
if ~isempty(MinPos)
    WordList(MinPos)=[];
    WordSize(MinPos)=[];
    WordNb(MinPos)=[];
end
if MaxSize>0
    MaxPos=find(WordSize>MaxSize);
    if ~isempty(MaxPos)
        WordList(MaxPos)=[];
        WordSize(MaxPos)=[];
        WordNb(MaxPos)=[];
    end
end

