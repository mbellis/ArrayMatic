%loads the last version of GEOmetadb
%=========================
% FUNCTION GEO_GETSQLITE %
%=========================
%
% GEO_GETSQLITE loads the last version of GEOmetadb from source at Meltzer' lab at NIH 
% modification of Yidong Chen, Genetics Branch, NCI, NIH, 2008 program
% 
%INPUT PARAMETERS
% 1 Action:
%       open geo metadb
%           open GEOmetadb.sqlite located in K.dir.geoMetadata
%       close
%           close GEO metadb
%       view tables and fields
%           view tables and fields existing in GEOmetadb.sqlite
%       write exemplar records
%       view technology
%       view species
%       display a GPL
%       print a GPL
%       import GSE
%       modif version
%       find biological conditions
%       import biological conditions
%
% FUNCTIONS
%
% VERIF
%
% PRINT
%
% EXTERNAL TOOL
% uses mksqlite.m developped by by Martin Kortmann <mail@kortmann.de>
% (sources downloaded at developer.berlios.de/projects/mksqlite/)

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

function geo_getsqlite()

global K
LogFileName = 'GEOmetadb.log';

% get current GEOmetadb time stamp.
WebPage = urlread( 'http://gbnci.abcc.ncifcrf.gov/geo/');
Tmp = regexp( WebPage, 'last_updated_GEOmetadb=(.*?) ', 'tokens' );
CurrTimeStamp = Tmp{1}{1};

% check whether there is an older version exists
cd(K.dir.geoMetadata)
if( exist( 'GEOmetadb.sqlite', 'file' ) )
    % if the db already exists, please check to see if it needs to be
    % updated.   
    % compare to timeStamp stored along with GEOmetadb.sqlite.
    if( exist( LogFileName, 'file' ) )
        [ImportDate, StoredTimeStamp] = textread( LogFileName, '%s%s', 'delimiter', '\t' );
        ImportDate=ImportDate{1};
        StoredTimeStamp=StoredTimeStamp{1};
        if( strcmp( StoredTimeStamp, CurrTimeStamp ) )
            % they are the same copy, exit!!
            h=msgbox(sprintf( 'No need to download (TimeStamp = %s). You have the latest GEOmetadb.sqlite.\n', CurrTimeStamp ));
            waitfor(h)
            return;
        end
    end
    
    % needs to be backup and download.
    movefile('GEOmetadb.sqlite',sprintf('GEOmetadb_%s.sqlite',ImportDate));
    fprintf( 'The old version of ''GEOmetadb.sqlite'' has been renamed to ''GEOmetadb_%s.sqlite''. \n',ImportDate );
end

%% find where this file sits:
% source of SQLite database file.
url = 'http://gbnci.abcc.ncifcrf.gov/geo/GEOmetadb.sqlite.gz';
gunzip(url,K.dir.geoMetadata);

DirInfo = dir('GEOmetadb.sqlite');
fprintf('The SQLite DB ''%s'' (size of %.1f MB) has been installed. \n\n', ...
    DirInfo.name, (DirInfo.bytes)/(1024*1024));

% generate log-file.
fid = fopen( fullfile( K.dir.geoMetadata, LogFileName ), 'w' );
fprintf( fid, '%s\t%s\n', date, CurrTimeStamp );
fclose( fid );
