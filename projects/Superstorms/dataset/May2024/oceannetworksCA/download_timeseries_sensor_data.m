function download_timeseries_sensor_data(varargin)
% function download_timeseries_sensor_data(my_token)
%    
% You will need an ONC personal token, which you can get from your user profile:
% from the 'Web Services API' tab at https://data.oceannetworks.ca/Profile when logged in.
%
% You can adjust the code in "call_api" to specify file "extension".  This
% code requests 'mat' files (as-is), but you can request csv, nc, or other
% formats. See the wiki for more information:
%   https://wiki.oceannetworks.ca/display/DP/1
%
% It can take a long time to prepare each data set, so you can create an
% extra filter called: partial_run = 1, which will quit after initiating
% the download. Then can either:
% 1) go and check your User Directory to manually download the mat files
% (click the "More" tab on the Oceans website, and select "User Directory")   
% 2) Call the API using the dpRunId that's been saved by the initial run
%
% The files you download contain metadata that includes a unique DOI for
% citation and reproducability.
%
% cd('F:\Documents\Projects\ADCP\ADCP Monitoring Management\SolarFlare\extended')

if nargin<1
    my_token = 'fbc71ba3-4984-4feb-a637-6e0af10524d8';  % Found in your ONC profile
else
    my_token = varargin{1};
end


%Download a month of data, around the May event
filters.dateFrom ='2024-05-01T00:00:01.000Z';
filters.dateTo = '2024-05-30T00:00:00.000Z';

%Can adjust this to select just a few data sets for download:
sites_ind = 1:18;
% sites_ind = [6,7,9,17,18];
% sites_ind = 1;

AA_codes =  {'ADCP1200KHZ';'ADCP150KHZ';'ADCP300KHZ';'ADCP600KHZ';'ADCP75KHZ';...
    'ADCP400KHZ';'ADCP55KHZ';'ADCP1MHZ';'ADCP2MHZ'};


for ii=sites_ind
    
    %try
    site = [];
    AA_code = [];
    
    if ii==1
        locationCode = 'BACAX';
        AA_code{1} = 'ADCP2MHZ';
        AA_code{2} = 'ADCP55KHZ';
        %24420 % Nortek Signature55 Current Profiler 200066
        %11302 % Nortek Aquadopp HR-Profiler 2 MHz 2700
    elseif ii==2
        locationCode = 'BACME';
        AA_code{1} = 'ADCP150KHZ';
        AA_code{2} = 'ADCP2MHZ';
        %11206 % RDI Workhorse Quartermaster ADCP 150 kHz (9455)
        %12003 % Nortek Aquadopp HR-Profiler 2 MHz 2978
    elseif ii==3
        locationCode = 'BACUS';
        AA_code{1} = 'ADCP600KHZ';
        AA_code{2} = 'ADCP2MHZ';
        %31399 % RDI Workhorse Monitor ADCP 600 kHz (SN 25472)
        %11203 % Nortek Aquadopp HR-Profiler 2 MHz ASP2699 AQD2973
    elseif ii==4
        locationCode = 'BIIP';
        AA_code{1} = 'ADCP600KHZ';
        %30699 % RDI Workhorse Monitor ADCP 600 kHz (SN 25471)
    elseif ii==5
        locationCode = 'BSM.J5';
        AA_code{1} = 'ADCP600KHZ';
        %73020 % RDI Workhorse Sentinel ADCP 600 kHz (SN 24287)
    elseif ii==6
        locationCode = 'FGPD';
        AA_code{1} = 'ADCP300KHZ';
        %20003 % RDI Workhorse Monitor ADCP 300 kHz (9273)
    elseif ii==7
        locationCode = 'FGPPN';
        AA_code{1} = 'ADCP600KHZ';
        AA_code{2} = 'ADCP2MHZ';
        %23498 % Nortek Aquadopp Profiler 2MHz 7401
        %23093 % RDI Workhorse Monitor ADCP 600 kHz (17434)
    elseif ii==8
        locationCode = 'HRBIP';
        AA_code{1} = 'ADCP2MHZ';
        %47140 % Nortek AWAC AST 400 kHz WAV 7749
    elseif ii==9
        locationCode = 'KEMF';
        AA_code{1} = 'ADCP600KHZ';
        %23063 % RDI Workhorse Monitor ADCP 600 kHz (17574)
    elseif ii==10
        locationCode = 'NC27';
        AA_code{1} = 'ADCP75KHZ';
        %23006 % RDI Workhorse Long Ranger ADCP 75 kHz (16511)
    elseif ii==11
        locationCode = 'NC89';
        AA_code{1} = 'ADCP75KHZ';
        %12108 % RDI Workhorse Long Ranger ADCP 75 kHz (9202)
    elseif ii==12
        locationCode = 'NCBC';
        AA_code{1} = 'ADCP75KHZ';
        %23065 % RDI Workhorse Long Ranger ADCP 75 kHz (17431)
    elseif ii==13
        locationCode = 'RCNE5';
        AA_code{1} = 'ADCP75KHZ';
        %13001 % RDI Workhorse Long Ranger ADCP 75 kHz (13429)
    elseif ii==14
        locationCode = 'RCNW5';
        AA_code{1} = 'ADCP75KHZ';
        %14001 % RDI Workhorse Long Ranger ADCP 75 kHz (3867)
    elseif ii==15
        locationCode = 'RCSE5.A1';
        AA_code{1} = 'ADCP75KHZ';
        %62780 % RDI Workhorse Long Ranger ADCP 75 kHz (24605)
    elseif ii==16
        locationCode = 'RCSE5.A2';
        AA_code{1} = 'ADCP75KHZ';
        %62760 % RDI Workhorse Long Ranger ADCP 75 kHz (24604)
    elseif ii==17
        locationCode = 'SCVIP';
        AA_code{1} = 'ADCP150KHZ';
        %23097 % RDI Workhorse Quartermaster ADCP 150 kHz (SN 17457)
    elseif ii==18
        locationCode = 'SEVIP';
        AA_code{1} = 'ADCP150KHZ';
        %65 % RDI Workhorse Quartermaster ADCP 150 kHz (SN 8580)
    end
    
    for jj=1:length(AA_code)
        filters.deviceCategoryCode = AA_code{jj};
        filters.locationCode = locationCode;
        % locationCode = 'BIIP';
        % deviceCategoryCode = 'ADCP600KHZ';
        filters.dataProductCode='TSSD';
        
        %Optional: Can download data manually through "User Directory"
        % filters.partial_run = 1;
        
        call_api(my_token, filters);
    end
    
    
end

end





function varargout = call_api(my_token, varargin)
% function call_onc_api(my_token, varargin)
% This is a bare-bones example created from the website:
%     https://wiki.oceannetworks.ca/display/DAQ/Request+Data+Product 
% As well as some code from onc\+onc\OncDelivery.m for downloading
% mat-files.
%
% This is probably the simplest way to download data through the ONC API
% using URLs 
%
% If requesting a type of device that is not set up yet, go to the DP
% documentation to find the correct data product code:
%  https://wiki.oceannetworks.ca/display/DP
% and search for text 'dataproductcode'
%
% my_token = '';  % Found in your ONC profile
% filters.locationCode = 'BIIP';
% filters.deviceCategoryCode = 'ADCP600KHZ';
% filters.dateFrom ='2024-03-01T00:00:00.000Z';
% filters.dateTo = '2024-03-02T01:00:00.000Z';
% call_onc_api(my_token, filters)
%
% Other possible filters: dataProductCode, average_interval (in seconds)
%
% Optional input: partial_run
% 1: run the data product request
% 2: run the data product download 
% backup option: use download_onc_ftp.m to download files from the ftp
% Ex: 
% Part1: Run the request. Will save a mat-file with the "run" structure 
% filters.partial_run = 1;
% run = call_onc_api(my_token, filters);
%
% Part2: Wait a bit, then go back and request the download, will load in a mat file with the "run" structure
% filters.partial_run = 2;
% call_onc_api(my_token, locationCode, filters);


%Initialize output:
varargout{1} = [];

%Define where data will be saved 
outPath = pwd;

%initialize parameters
locationCode=[]; deviceCategoryCode=[];  dateFrom=[];  dateTo=[];  dataProductCode=[]; average_interval = []; partial_run = 0;

if nargin >1
    filters = varargin{1};
    locationCode = filters.locationCode;
    deviceCategoryCode = filters.deviceCategoryCode;
    dateFrom = filters.dateFrom;
    dateTo = filters.dateTo;
    
    %Check if a specific data product code is requested, otherwise will
    %choose a default based on the device
    if isfield(filters, 'dataProductCode')
        if ~isempty(filters.dataProductCode)
            dataProductCode = filters.dataProductCode;
        end
    end
    
    %Check if averaging is specified, otherwise defaults to no averaging
    if isfield(filters, 'average_interval')
        if ~isempty(filters.average_interval)
            average_interval = filters.average_interval;
        end
    end
    
    %Check if a partial run was requested:
    if isfield(filters, 'partial_run')
        if ~isempty(filters.partial_run)
            partial_run = filters.partial_run;
        end
    end
    
end

    

if partial_run==0 || partial_run==1
    %***************************************************
    %PART 1: returns a Request Id, which can be used to run the data product.
    %***************************************************
    requestInfo = initiate_request_id(my_token, locationCode, deviceCategoryCode, dateFrom, dateTo, dataProductCode, average_interval);
    
    %***************************************************
    %PART 2: runs a data product request using a Request Id returned from the
    %'Request Data Product' example
    %***************************************************
    run = initiate_data_request(my_token, requestInfo);
    
    %If desired, return the "run" information and quit this code
    if partial_run == 1
        varargout{1} = run;
        dateFrom_str = datestr(datenum(dateFrom,'yyyy-mm-ddTHH:MM:SS.FFFZ'),'yyyymmddTHHMMSS');
        dateTo_str = datestr(datenum(dateTo,'yyyy-mm-ddTHH:MM:SS.FFFZ'),'yyyymmddTHHMMSS');
        run_info_file = fullfile(pwd,sprintf('%s_%s_%s_%s_runinfo.mat',locationCode,deviceCategoryCode, dateFrom_str, dateTo_str));
        save(run_info_file,'run')
        return;
    else
        if nargout == 1
            varargout{1} = [];
        end
    end
    
end

if partial_run==0 || partial_run==2
    %***************************************************
    %PART 3: Download all of the files generated by a data product request
    %using a Run Id returned from the 'Run Data Product Request' example
    %***************************************************
    
    %Load the "run" structure from file
    if partial_run == 2
        dateFrom_str = datestr(datenum(dateFrom,'yyyy-mm-ddTHH:MM:SS.FFFZ'),'yyyymmddTHHMMSS');
        dateTo_str = datestr(datenum(dateTo,'yyyy-mm-ddTHH:MM:SS.FFFZ'),'yyyymmddTHHMMSS');
        run_info_file = fullfile(pwd,sprintf('%s_%s_%s_%s_runinfo.mat',locationCode,deviceCategoryCode, dateFrom_str, dateTo_str));
        if exist(run_info_file, 'file')
            load(run_info_file);
        else
           error('could not find file named : %s \n',run_info_file); 
        end
    end
        
    download_data_product(my_token, run, outPath)
end




end



function requestInfo = initiate_request_id(my_token, locationCode, deviceCategoryCode, dateFrom, dateTo, dataProductCode, average_interval)
    %***************************************************
    %PART 1: returns a Request Id, which can be used to run the data product.
    %***************************************************

    if nargin == 1
        %Default example from the API website:
        % https://wiki.oceannetworks.ca/display/O2A/Request+Data+Product
        url = ['https://data.oceannetworks.ca/api/dataProductDelivery', ...
            '?method=request' ...
            '&token=',my_token ...                            % replace YOUR_TOKEN_HERE with your personal token obtained from the 'Web Services API' tab at https://data.oceannetworks.ca/Profile when logged in.
            '&locationCode=BACAX' ...                               % Barkley Canyon / Axis (POD 1)
            '&deviceCategoryCode=ADCP2MHZ' ...                      % 150 kHz Acoustic Doppler Current Profiler
            '&dataProductCode=TSSD' ...                             % Time Series Scalar Data
            '&extension=csv' ...                                    % Comma Separated spreadsheet file
            '&dateFrom=2016-07-27T00:00:00.000Z' ...                   % The datetime of the first data point (From Date)
            '&dateTo=2016-08-01T00:00:00.000Z' ...                     % The datetime of the last data point (To Date)
            '&dpo_qualityControl=1' ...                             % The Quality Control data product option - See https://wiki.oceannetworks.ca/display/DP/1
            '&dpo_resample=none' ...                                % The Resampling data product option - See https://wiki.oceannetworks.ca/display/DP/1
            '&dpo_dataGaps=0'];                                     % The Data Gaps data product option - See https://wiki.oceannetworks.ca/display/DP/1

        request = matlab.net.http.RequestMessage;
        uri = matlab.net.URI(url);

    elseif strcmpi(dataProductCode,'TSSD')

        %Define the base url
        base_url = 'https://data.oceannetworks.ca/api/dataProductDelivery';
        extension = 'mat'; %Options include: rdi, mat, nc
        %Optional: Use user-specified averaging: Other options [s]: 60, 300, 600, 900, 3600
        dpo_average = 60;
        if ~isempty(average_interval)
            dpo_average = average_interval;
        end
        %dpo_average = sprintf('%d',dpo_average); %Not sure this is needed
        %ensemblePeriod = 0; %Default is no averaging.


        filters_ = {'method','request' ...
            'token',my_token,...
            'locationCode', locationCode,...
            'deviceCategoryCode', deviceCategoryCode,...
            'dataProductCode', dataProductCode,...
            'extension', extension,...
            'dateFrom', dateFrom,...
            'dateTo', dateTo,...                     % The datetime of the last data point (To Date)
            'dpo_qualityControl', '1',...                             % The Quality Control data product option - See https://wiki.oceannetworks.ca/display/DP/1
            'dpo_resample','average',...%'none' ...                                % The Resampling data product option - See https://wiki.oceannetworks.ca/display/DP/1
            'dpo_average',dpo_average,... %This is only included if dpo_resample is 'average'
            'dpo_dataGaps','1'};                                     % The Data Gaps data product option - See https://wiki.oceannetworks.ca/display/DP/1

        %             'dpo_qualityControl', '0'...                             % The Quality Control data product option - See https://wiki.oceannetworks.ca/display/DP/1
        %             'dpo_resample','none'...
        %             'dpo_dataGaps','1'};                                     % The Data Gaps data product option - See https://wiki.oceannetworks.ca/display/DP/1

        % sanitize filters
        filters = sanitize_filters(filters_);

        % prepare HTTP request
        request = matlab.net.http.RequestMessage;
        uri = matlab.net.URI(base_url);
        uri.Query = matlab.net.QueryParameter(filters);
        %url = char(uri);


    else
        %Define the base url
        base_url = 'https://data.oceannetworks.ca/api/dataProductDelivery';
        extension = 'mat'; %Options include: rdi, mat, nc
        %Optional: Use user-specified averaging: Other options [s]: 60, 300, 600, 900, 3600
        ensemblePeriod = 0; %Default is no averaging.
        if ~isempty(average_interval)
            ensemblePeriod = average_interval;
        end
        %dpo_average = sprintf('%d',dpo_average); %Not sure this is needed

        if ismember(deviceCategoryCode,{'ADCP75KHZ','ADCP150KHZ','ADCP300KHZ','ADCP600KHZ'})
            %Assume we want ADCP data
            dataProductCode = 'RADCPTS'; %'RDI ADCP Time Series'

            %See the help doc for details about all filters:
            %helpDocument: 'https://wiki.oceannetworks.ca/display/DP/5'

            filters_ = {'method','request' ...
                'token',my_token,...
                'locationCode', locationCode,...
                'deviceCategoryCode', deviceCategoryCode,...
                'dataProductCode', dataProductCode,...
                'extension', extension,...
                'dateFrom', dateFrom,...
                'dateTo', dateTo,...
                'dpo_ensemblePeriod', ensemblePeriod,... % 0 = data not altered (none). All others in s: 60, 600, 900, 3600
                'dpo_velocityBinmapping',  -1,... % -1 = As configured on device.  1 = Nearest, 0 = None, 2 = Linear interpolation
                'dpo_3beam', 'config',... %config (as configured on device), On,  Off
                'dpo_corScreen', -1,... %-1 = as configured on device,  All other in counts: 64, 0, 16, 32, 128
                'dpo_errVelScreen', -1,... %-1 = as configured on device.  All others in m/s : 2, 0, 5, 1, 0.5, 0.25, 0.1
                'dpo_falseTarScreen', -1,... %-1 = as configured on device.  255, 192, 128, 64, 50, 32, 16
                %'deleteFile','false'
                };

            %Optional if want 3-beam on:
            %  if do_3beam
            %     filters_(18) = {'On'};%config (as configured on device), On,  Off
            %  end

            do_test = 0;
            if do_test
                %Test for Rich Pawlowich
                filters_ = {'method','request' ...
                    'token',my_token,...
                    'locationCode', locationCode,...
                    'deviceCategoryCode', deviceCategoryCode,...
                    'dataProductCode', dataProductCode,...
                    'extension', extension,...
                    'dateFrom', '2024-03-17T00:00:00.000Z',...
                    'dateTo', '2024-03-17T23:59:59.999Z',...
                    'dpo_ensemblePeriod', ensemblePeriod,...
                    'dpo_velocityBinmapping', 1,...
                    'dpo_3beam', 'On',...
                    'dpo_corScreen', 0,...
                    'dpo_errVelScreen', 0,...
                    'dpo_falseTarScreen', 255,...
                    };
            end

        elseif ismember(deviceCategoryCode,{'ADCP400KHZ','ADCP1MHZ','ADCP2MHZ'})
            dataProductCode = 'NTS'; %'Nortek Time Series'

            %See the help doc for details about all filters:
            %helpDocument: 'https://wiki.oceannetworks.ca/display/DP/22'

            filters_ = {'method','request' ...
                'token',my_token,...
                'locationCode', locationCode,...
                'deviceCategoryCode', deviceCategoryCode,...
                'dataProductCode', dataProductCode,...
                'extension', extension,...
                'dateFrom', dateFrom,...
                'dateTo', dateTo,...
                'dpo_ensemblePeriod', ensemblePeriod,... % 0 = data not altered (none). All others in s: 60, 600, 900, 3600
                'dpo_velocityBinmapping',  1,... % -1 = As configured on device.  1 = Nearest, 0 = None, 2 = Linear interpolation
                %'deleteFile','false'
                %'corScreen', 0,... %50 is the default,  All other in percents: 0,5,10, 20, 40, 60, 70, 90, 95 (NOTE: Invalid for AWAC)
                %'dpo_nortekSignatureSeriesPlan',0 %For nortek55's only! 0 or 1
                };

        elseif ismember(deviceCategoryCode,{'ADCP55KHZ'})
            dataProductCode = 'NTS'; %'Nortek Time Series'

            %See the help doc for details about all filters:
            %helpDocument: 'https://wiki.oceannetworks.ca/display/DP/22'

            filters_ = {'method','request' ...
                'token',my_token,...
                'locationCode', locationCode,...
                'deviceCategoryCode', deviceCategoryCode,...
                'dataProductCode', dataProductCode,...
                'extension', extension,...
                'dateFrom', dateFrom,...
                'dateTo', dateTo,...
                'dpo_ensemblePeriod', ensemblePeriod,... % 0 = data not altered (none). All others in s: 60, 600, 900, 3600
                'dpo_velocityBinmapping',  1,... % -1 = As configured on device.  1 = Nearest, 0 = None, 2 = Linear interpolation
                'dpo_nortekSignatureSeriesPlan', 'both',... %For nortek55's only! 0 or 1, or 'all' for both
                'dpo_corScreen', 0,... %50 is the default,  All other in percents: 0,5,10, 20, 40, 60, 70, 90, 95 (NOTE: Invalid for AWAC)
                %'deleteFile','false'
                };

        elseif ismember(deviceCategoryCode,{'CURRENTMETER'})
            dataProductCode = 'NTS'; %'Nortek Time Series'

            %special case: don't get the mat file:
            extension = 'aqd'; %Assumes this is raw format (vec, prf, aqd, etc)

            %See the help doc for details about all filters:
            %helpDocument: 'https://wiki.oceannetworks.ca/display/DP/22'

            filters_ = {'method','request' ...
                'token',my_token,...
                'locationCode', locationCode,...
                'deviceCategoryCode', deviceCategoryCode,...
                'dataProductCode', dataProductCode,...
                'extension', extension,...
                'dateFrom', dateFrom,...
                'dateTo', dateTo,...
                %'dpo_ensemblePeriod', ensemblePeriod,... % 0 = data not altered (none). All others in s: 60, 600, 900, 3600
                %'dpo_velocityBinmapping',  1,... % -1 = As configured on device.  1 = Nearest, 0 = None, 2 = Linear interpolation
                %'dpo_nortekSignatureSeriesPlan', 'both',... %For nortek55's only! 0 or 1, or 'all' for both
                %'dpo_corScreen', 0,... %50 is the default,  All other in percents: 0,5,10, 20, 40, 60, 70, 90, 95 (NOTE: Invalid for AWAC)
                %'deleteFile','false'
                };
            
        elseif ismember(deviceCategoryCode,{'ECHOSOUNDERBIOA'})
            dataProductCode = 'AAPTS'; %'Nortek Time Series'
            
            %special case: don't get the mat file:
            extension = '01a'; 
            %extension = 'xml'; 
            
            %See the help doc for details about all filters:
            %helpDocument: 'https://wiki.oceannetworks.ca/display/DP/24'
            
            filters_ = {'method','request' ...
                'token',my_token,...
                'locationCode', locationCode,...
                'deviceCategoryCode', deviceCategoryCode,...
                'dataProductCode', dataProductCode,...
                'extension', extension,...
                'dateFrom', dateFrom,...
                'dateTo', dateTo,...
                %'dpo_ensemblePeriod', ensemblePeriod,... % 0 = data not altered (none). All others in s: 60, 600, 900, 3600
                %'dpo_velocityBinmapping',  1,... % -1 = As configured on device.  1 = Nearest, 0 = None, 2 = Linear interpolation
                %'dpo_nortekSignatureSeriesPlan', 'both',... %For nortek55's only! 0 or 1, or 'all' for both
                %'dpo_corScreen', 0,... %50 is the default,  All other in percents: 0,5,10, 20, 40, 60, 70, 90, 95 (NOTE: Invalid for AWAC)
                %'deleteFile','false'
                };

        end


        % sanitize filters
        %I got this function from: Ocean Networks Canada API Client Library\onc\+util\sanitize_filters.m
        %It takes the matlab struct and converts it into a string
        filters = sanitize_filters(filters_);

        % prepare HTTP request
        request = matlab.net.http.RequestMessage;
        uri = matlab.net.URI(base_url);
        uri.Query = matlab.net.QueryParameter(filters);
        %url = char(uri);

    end
    %Submit the request, and get info on how large it is, and get a requestId
    options = matlab.net.http.HTTPOptions('ConnectTimeout',60);
    response = send(request,uri,options);

    %Check the statusCode of the response
    if (response.StatusCode == 200)    % HTTP Status - OK
        requestInfo = response.Body.Data;
        disp(requestInfo)
    elseif (response.StatusCode == 400) % HTTP Status - Bad Request
        error(response.Body.Data.errors);
    else % all other HTTP Statuses
        disp(char(response.StatusLine));
    end
    
    %Add extension to the requestInfo
    requestInfo.extension = extension;

end


function run = initiate_data_request(my_token, requestInfo)
    %***************************************************
    %PART 2: runs a data product request using a Request Id returned from the
    %'Request Data Product' example
    %***************************************************
    url = ['https://data.oceannetworks.ca/api/dataProductDelivery', ...
        '?method=run' ...
        '&token=',my_token ...                    % replace YOUR_TOKEN_HERE with your personal token obtained from the 'Web Services API' tab at https://data.oceannetworks.ca/Profile  when logged in.
        '&dpRequestId=',num2str(requestInfo.dpRequestId)];           % replace YOUR_REQUEST_ID_HERE with a requestId number returned from the request method

    request = matlab.net.http.RequestMessage;
    uri = matlab.net.URI(url);
    options = matlab.net.http.HTTPOptions('ConnectTimeout',60);

    response = send(request,uri,options);
    % statusLine = split(response.StatusLine,' ')
    % statusCode = response.StatusLine(


    if any(response.StatusCode == [200,202])  % HTTP Status is 200 - OK or 202 - Accepted
        runs = response.Body.Data;
        for i=1:numel(runs)
            run = runs(i);
            if (isfield(run,'queuePosition'))
                disp(sprintf('Run Id: %i, Status: %s, Queue Position: %s',run.dpRunId,run.status,num2str(run.queuePosition)));
            else
                disp(sprintf('Run Id: %i, Status: %s',run.dpRunId,run.status));
            end
        end
        %Add extension to the run structure
        run.extension = requestInfo.extension;
        
    elseif (response.StatusCode == 400) % HTTP Status - Bad Request
        errors = response.Body.Data.errors;
        for i=1:numel(errors)
            disp(errors(i));
        end
    else % all other HTTP Statuses
        disp(char(response.StatusLine));
    end
end


function download_data_product(my_token, run, outPath)
    %***************************************************
    %PART 3: Download all of the files generated by a data product request
    %using a Run Id returned from the 'Run Data Product Request' example
    %***************************************************
    index_num = 1;  %Specify first file for download
    [uri, response] = create_url(my_token, run, index_num);

    %IDEA: It might be a good idea to put all of the code below into a seperate
    %function, and put that inside of the "while doLoop".

    %Check if a message about the progress is available. Continue to check the
    %messages until something other than 202 occurs.
    message = ''; pollPeriod = 1;
    while response.StatusCode == 202
        payload = response.Body.Data;
        if (isfield(payload,'message'))
            message_new = sprintf('HTTP %i - %s: %s',response.StatusCode,getReasonPhrase(response.StatusCode),payload.message');
            if ~strcmpi(message, message_new)
                fprintf('\n');
                disp(message_new);
                message = message_new;
            else
                fprintf('.');
            end
        else
            disp(payload);
        end

        pause(pollPeriod); %pause for 1 second;

        %Call again, see if it's been updated
        [uri, response] = create_url(my_token, run, index_num);
    end
    fprintf('\n')



    %Now that the statusCode is no longer 202, Data is ready to download:
    if (response.StatusCode == 200)        % HTTP OK

        if strcmpi(run.extension,'mat') || strcmpi(run.extension,'01a') || strcmpi(run.extension,'xml')
            %I tested this on .01a and .xml data (for AZFPs) and it worked
            %just fine.  I'm not sure why a seperate method exists. 
            
            % run.dpRunId = 38994936; %For testing

            % file was downloaded (200), or downloaded & skipped (777)

            doLoop = true;
            while doLoop
                %This is from onc\+onc\DataProductFile, L 148
                if response.StatusCode == 200

                    %Parse the filename
                    txt = response.getFields('Content-Disposition').Value;
                    tokens = split(txt, 'filename=');
                    fileName = tokens(2);
                    %filename = this.extractNameFromHeader(response);
                    [~, saveResult] = urlwrite(uri.EncodedURI,fullfile(outPath, fileName));

                    % log status
                    if saveResult == 0
                        fprintf('Downloaded "%s" \n',fileName);
                    elseif saveResult == -2
                        fprintf('Skipping "%s": File already exists\n', fileName);
                        endStatus = 777;
                    end

                    %I haven't tested any of these! Don't really know how
                    %they'll be handled. I've commented out anything that I
                    %think will cause issues.
                elseif response.StatusCode == 202
                    % Still processing, wait and retry
                    %log.printResponse(jsondecode(response.Body.Data));
                    pause(pollPeriod);
                elseif response.StatusCode == 204
                    % No data found
                    fprintf('No data found.\n');
                elseif response.StatusCode == 400
                    % API Error
                    %util.print_error(response, fullUrl);
                    err = struct( ...
                        'message', sprintf('The request failed with HTTP status %d\n', response.StatusCode), ...
                        'identifier', 'DataProductFile:HTTP400');
                    error(err);
                elseif response.StatusCode == 404
                    % Index too high, no more files to download
                else
                    % File is gone
                    fprintf('ERROR: File with runId %d and index "%s" not found\n', ...
                        run.dpRunId, index_num);
                    %util.print_error(response, fullUrl);
                end


                %Increment the index and create a new download uri, which
                %will be used on the next loop to download the next file.
                index_num = index_num+1;
                [uri, response] = create_url(my_token, run, index_num);

                %If 'response' indicates that no more files exist, exit the loop
                if response.StatusCode~=200 %endStatus ~= 202
                    % no more files to download
                    doLoop = false;
                end

            end


        else
            %DEFAULT: If downloading csv or other data:
            %This fails for xml and mat file types.  
            doLoop = true;
            while doLoop
                if response.StatusCode == 200
                    fileName = '';
                    size = 0;
                    for i=1:numel(response.Header)
                        fld = response.Header(i);
                        if (fld.Name == "Content-Disposition")
                            S = strsplit(fld.Value,'=');
                            fileName = S(2);
                        end
                        if (fld.Name == "Content-Length")
                            size = fld.Value;
                        end
                    end
                    
                    file = sprintf('%s\\%s',outPath,fileName);
                    fprintf('writing file: %s\n',file);
                    fileID = fopen(file,'w');
                    cleanID = onCleanup(@() fclose(fileID));
                    fwrite(fileID,response.Body.Data);
                end
                
                %Increment the index and create a new download uri, which
                %will be used on the next loop to download the next file.
                index_num = index_num+1;
                [uri, response] = create_url(my_token, run, index_num);

                %If 'response' indicates that no more files exist, exit the loop
                if response.StatusCode~=200 %endStatus ~= 202
                    % no more files to download
                    doLoop = false;
                end
                
            end
        end

    elseif (any(response.StatusCode == [202,204,404,410]))
        payload = response.Body.Data;
        if (isfield(payload,'message'))
            disp(sprintf('HTTP %i - %s: %s',response.StatusCode,getReasonPhrase(response.StatusCode),payload.message'));
        else
            disp(payload);
        end
    elseif (r.StatusCode == 400) % HTTP BadRequest
        errors = response.Body.Data.errors;
        for i=1:numel(errors)
            disp(errors(i));
        end
    else % all other HTTP Statuses
        disp(response.Body.Data);
    end

end


function [uri, response] = create_url(my_token, run, index_num)
    %Create the url for downloading data. This is placed in a function
    %because it must be created until all files downloaded, and only way to
    %know if all files have been downloaded is if an error is returned
    %because the index_num is too large. 
    url = ['https://data.oceannetworks.ca/api/dataProductDelivery' ...
        '?method=download' ...
        '&token=',my_token ...                         % replace YOUR_TOKEN_HERE with your personal token obtained from the 'Web Services API' tab at https://data.oceannetworks.ca/Profile when logged in..
        '&dpRunId=',num2str(run.dpRunId) ...           % replace YOUR_RUN_ID with the dpRunId returned from the 'run' method.
        '&index=',num2str(index_num)];                 % for run requests that contain more than one file, change the index number to the index of the file you would like to download.
    % If the index number does not exist an HTTP 410 and a message will be returned.
    request = matlab.net.http.RequestMessage;
    uri = matlab.net.URI(url);
    options = matlab.net.http.HTTPOptions('ConnectTimeout',60);
    response = send(request,uri,options);
end


% Preprocess filters to fit the expected structure
function r = sanitize_filters(filters)
    % Make sure filters are a struct
    if class(filters) == "cell"
        r = struct();
        
        s = size(filters);
        rows = s(1);
        cols = s(2);
        
        if (rows > 1) && (cols == 2)
            % if the cell array has 2 columns, interpret properly
            for row = 1 : rows
                name = filters{row, 1};
                value = filters{row, 2};
                r.(name) = value;
            end
        else
            % suppose filters is a 1-row cell array
            n = numel(filters);
            for i = 1 : 2 : n
                name = filters{i};
                value = filters{i + 1};
                r.(name) = value;
            end
        end
    else
        r = filters;
    end
    
    % translate boolean values to strings
    names = fieldnames(r);
    for i = 1 : numel(names)
        name = names{i};
        if strcmp(class(r.(name)), "logical")
            if (r.(name) == true)
                r.(name) = "true";
            else
                r.(name) = "false";
            end
        end
    end
end




