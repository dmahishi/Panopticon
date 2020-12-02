function [h,v] = PanoPlot(objlist,paramin)
%PANOPLOT Ploting function for panopticon experiments.
%   Use PanoPlot to create plots
%   PanoPlot(objlist,paramin) creates plots for the datasets in objlist (a
%   cell array). paramin contains the relevant parameters (also cell array)
%
%   Examples:
%   PanoPlot({pdhlA,pdhlB},{'defaultPlot','Activity', 'binSize', 1800}
%
%   Valid parameters:
%   defaultPlot: use default settings for standardized plots
%       'Activity','PatchActivity','ActivityGroups',
%       'StopDistributionBasic','StopDistributionMulti','MovementDistribution'
%   binSize: used for setting the size of time bins 
%       default: 1800
%   xTicks: to define x-axis tick number series
%       default: 2
%   subPlotRange: Select the max/min time range for trimming and plotting
%       default: []
%   legendLabel: To assign field-specific legend labels from Metadata info
%       default: Substrate (eg. Sucrose, Peptone, Agarose etc..)
%   plotParams:  which specific area of the Arena to plot data from
%       default: (1) {(2) for food patches}
%
%   Authors: Tilman Triphan, tilman.triphan@uni-leipzig.de
%            Deepthi Mahishi, deepthi.mahishi_vasuki@uni-leipzig.de
%   License: GNU General Public License v3.0
%   Copyright (C) 2020  Tilman Triphan, Deepthi Mahishi
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   This functions requires raacampbell/shadedE​rrorBar to work.
%   www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar

p = inputParser;

isDefaultPlot = @(x) any(validatestring(x,{'Activity','PatchActivity',...
    'ActivityGroups','StopDistributionBasic','StopDistributionMulti',...
    'MovementDistribution'}));

isPlotFunc = @(x) any(validatestring(x,{'isActive','isActiveStacked',...
    'isActiveStacked','onArea','onAreaAll','farFromArea',...
    'nearArea','closeToOrOnArea','distance', 'distanceToArea'}));

isPlotType = @(x) any(validatestring(x,{'shadedErrorBar','boxPlots',...
    'stackedBars', 'bar'}));

isBinMethod = @(x) any(validatestring(x,{'auto','mean','sum','log2'}));

isValidBin = @(x) ismember(x,[30 60 120 300 600 900 1800 3600 7200 10800 14400]);

isInfoField =  @(x) any(validatestring(x,fieldnames(objlist{1}(1).getInfo)));

defaultColors = {[0.00 0.45 0.74],[0.05 0.37 0.07],[0.64 0.08 0.18]};

p.KeepUnmatched = 1;
p.addParameter('defaultPlot','',isDefaultPlot);
p.parse(paramin{:});
params = p.Results;

switch (params.defaultPlot)
    % Please include metadata infofield exact names
    % Please refer above for the list of allowed params for each case.
    case 'Activity'
        p.addParameter('toPlot','isActive',isPlotFunc); 
        p.addParameter('plotType','shadedErrorBar',isPlotType);
        params.indices = ':';
        p.addParameter('binSize', 1800,isValidBin);
        p.addParameter('binMethod','mean',isBinMethod);
        p.addParameter('multiRounds',false,@(x) islogical(x));
        p.addParameter('xTicks',2,@(x) isnumeric(x) && (x > 0));
        p.addParameter('subPlotRange',[],@(x) isnumeric(x) && all(x > 0));
        p.addParameter('legendLabel','Starvation',isInfoField);
        p.addParameter('plotParams',[],@(x) iscell(x));
        p.addParameter('saveDataAs','',@(x) ischar(x)); % works only for single datasets        
    case 'PatchActivity'
        p.addParameter('toPlot','onArea',isPlotFunc);
        p.addParameter('plotType','shadedErrorBar',isPlotType);
        p.addParameter('indices',':', @(x) isnumeric(x) || islogical(x));
        p.addParameter('binSize',1800,isValidBin);
        p.addParameter('binMethod','mean',isBinMethod);
        p.addParameter('multiRounds',false,@(x) islogical(x));
        p.addParameter('xTicks',2,@(x) isnumeric(x) && (x > 0));
        p.addParameter('subPlotRange',[],@(x) isnumeric(x) && all(x > 0));
        p.addParameter('legendLabel','Starvation',isInfoField);
        p.addParameter('plotParams',{2},@(x) iscell(x));
        p.addParameter('saveDataAs','',@(x) ischar(x)); % works only for single datasets 
    case 'ActivityGroups'
        p.addParameter('toPlot','isActiveStackedMM',isPlotFunc);
        p.addParameter('plotType','stackedBars',isPlotType);
        p.addParameter('indices',':',@(x) isnumeric(x) || islogical(x));
        p.addParameter('binSize',1800,isValidBin);
        p.addParameter('binMethod','mean',isBinMethod);
        p.addParameter('multiRounds',true,@(x) islogical(x)); 
        p.addParameter('xTicks',2,@(x) isnumeric(x) && (x > 0));
        p.addParameter('subPlotRange',[],@(x) isnumeric(x) && all(x > 0));
        p.addParameter('plotParams',{'AT(a)','AT(a+1)'},@(x) iscell(x)); % Activity thresholds can be set in pdhl
    case 'StopDistributionBasic'
        p.addParameter('toPlot','isActive',isPlotFunc);
        p.addParameter('plotType','boxPlots',isPlotType);
        p.addParameter('indices',':',@(x) isnumeric(x) || islogical(x));
        p.addParameter('timeWindow',[],@(x) isnumeric(x) && all(x >= 0) && (numel(x) == 2)); %[8.25 12.25] for 8:15 to 12:15        
        p.addParameter('binCount',10,@(x) isnumeric(x) && (x > 0));
        p.addParameter('binWidth',10,@(x) isnumeric(x) && (x > 0));
        p.addParameter('binMethod','log2',isBinMethod);
        p.addParameter('includeHighValues',true,@(x) islogical(x));
        p.addParameter('keepMinimum',2,@(x) isnumeric(x) && (x > 0));
        p.addParameter('plotParams',{1},@(x) iscell(x));
        p.addParameter('saveDataAs','',@(x) ischar(x)); % works only for single datasets 
    case 'StopDistributionMulti'
        p.addParameter('toPlot','isActive',isPlotFunc);
        p.addParameter('plotType','boxPlots',isPlotType);
        p.addParameter('indices',':',@(x) isnumeric(x) || islogical(x));
        p.addParameter('timeWindow',[],@(x) isnumeric(x) && all(x >= 0) && (numel(x) == 2)); %[8.25 12.25] for 8:15 to 12:15        
        p.addParameter('binCount',10,@(x) isnumeric(x) && (x > 0));
        p.addParameter('binWidth',10,@(x) isnumeric(x) && (x > 0));
        p.addParameter('binMethod','log2',isBinMethod);
        p.addParameter('includeHighValues',true,@(x) islogical(x));
        p.addParameter('keepMinimum',2,@(x) isnumeric(x) && (x > 0));
        p.addParameter('legendLabel','Starvation',isInfoField);
        p.addParameter('plotParams',{1},@(x) iscell(x));
    case 'MovementDistribution' 
        p.addParameter('toPlot','distance',isPlotFunc);
        p.addParameter('plotType','histogram',isPlotType);                  % type of plot
        p.addParameter('indices',':',@(x) isnumeric(x) || islogical(x));    % range of data to plot
        p.addParameter('timeWindow',[],@(x) isnumeric(x) && all(x >= 0) && (numel(x) == 2)); %[8.25 12.25] for 8:15 to 12:15
        p.addParameter('binCount',75,@(x) isnumeric(x) && (x > 0));        % total number of bins 
        p.addParameter('binWidth',0.4,@(x) isnumeric(x) && (x > 0));        
        p.addParameter('binMethod','auto',isBinMethod);                     
        p.addParameter('binLimits',[0 15],@(x) isnumeric(x) && all(x >= 0) && (numel(x) == 2));
        p.addParameter('keepMinimum',0.2,@(x) isnumeric(x) && (x >= 0));    %
        p.addParameter('normalization','probability',@(x) ischar(x));       % check Matlab histogram function for valid normalization schemes
        p.addParameter('filterFunction',[],@(x) ischar(x));                 % movement distribution on food patch ('filterFunction', 'onArea')
        p.addParameter('filterParams',{1},@(x) iscell(x));                  % movement distribution on food patch ('filterParams', {2})
        p.addParameter('legendLabel','Starvation',isInfoField);
        p.addParameter('plotParams',{1},@(x) iscell(x));
    otherwise
        p.addParameter('toPlot','',isPlotFunc);
        p.addParameter('plotType','',isPlotType);
        p.addParameter('indices',':',@(x) isnumeric(x) || islogical(x));
        p.addParameter('binSize',1800,isValidBin);
        p.addParameter('binCount',10,@(x) isnumeric(x) && (x > 0));
        p.addParameter('binMethod','mean',isBinMethod);
        p.addParameter('includeHighValues',true,@(x) islogical(x));
        p.addParameter('keepMinimum',2,@(x) isnumeric(x) && (x > 0));
        p.addParameter('multiRounds',false,@(x) islogical(x));
        p.addParameter('xTicks',2,@(x) isnumeric(x) && (x > 0));
        p.addParameter('subPlotRange',[],@(x) isnumeric(x) && all(x > 0));
        p.addParameter('legendLabel','Substrate',isInfoField);
        p.addParameter('plotParams',[],@(x) iscell(x));
end

p.addParameter('fontSize',25,@(x) isnumeric(x) && (x > 0));
p.addParameter('colorList',defaultColors,@(x) isnumeric(x));

p.parse(paramin{:});
params = p.Results;

h = figure;
hold on

if strcmp(params.defaultPlot,'StopDistributionBasic')
    for a=1:numel(objlist)
        obj = objlist{a};
        bin = zeros(numel(obj),params.binCount);
        for o=1:numel(obj)
            [posEventStarts, durEvents] = obj(o).calculateStops(params.indices);

            lOnArea = obj(o).onArea(params.plotParams{1},posEventStarts);
            durEvents = durEvents(lOnArea);
            
            idx = durEvents < params.keepMinimum;
            durEvents(idx) = [];
            if ~isempty(durEvents)
                durEvents(length(durEvents)) = [];
            end
            for i=1:params.binCount
                switch (params.binMethod)
                    case 'log2'
                        bin(o,i) = sum((durEvents >= 2^i) & (durEvents < 2^(i+1)-1));
                        if (i == params.binCount) && params.includeHighValues
                            bin(o,i) = sum(durEvents >= 2^i);
                        end
                    otherwise
                        bin(o,i) = sum(sum(durEvents == (i-1)*params.binWidth+1:i*params.binWidth));
                        if (i == params.binCount) && params.includeHighValues
                            bin(o,i) = sum(durEvents >= (i-1)*params.binWidth+1);
                        end
                end
            end
        end
    end
    if ~isempty(params.saveDataAs)
        %Saves data along with exp ID labels
        cBinVals = arrayfun(@(x) x, bin, 'UniformOutput', false);
        xlsdata = cell(numel(obj), params.binCount+1);
        xlsdata(:,1) = obj.getInfo('SectorID');
        xlsdata(:,2:end) = cBinVals;
        xlswrite(params.saveDataAs,xlsdata);
    end
    boxplot(bin);
    axis([0.5 10.5 0 2000]);
    switch (params.binMethod)
        case 'log2'
            xlabel('Stop duration(log2)','FontSize',params.fontSize,'FontWeight','bold');
        otherwise
            xlabel(['Stop duration(',num2str(params.binWidth),'s/bin)'],...
                'FontSize',params.fontSize,'FontWeight','bold');
    end
    ylabel('Number of Stop Events','FontSize',params.fontSize,'FontWeight','bold')
    return
end

%--------------------------------------------------------------------------
if strcmp(params.defaultPlot,'StopDistributionMulti')
    for a=1:numel(objlist)
        obj = objlist{a};
        bin = zeros(numel(obj),params.binCount);
        for o=1:numel(obj)
            %if time indices, calculate here            
            if ~isempty(params.timeWindow)
                tws = obj(o).timeToFrame(params.timeWindow(1));
                twe = obj(o).timeToFrame(params.timeWindow(2));
                params.indices = tws:twe;
            end
            [posEventStarts, durEvents] = obj(o).calculateStops(params.indices);
            if posEventStarts > 0
                onArea = obj(o).onArea(params.plotParams{1},params.indices);
                stopOnAreaStatus = onArea(posEventStarts);
                durEventsOnArea = durEvents(stopOnAreaStatus);
                durEvents = durEventsOnArea;
            else
                durEvents = [];
            end
            
            idx = durEvents < params.keepMinimum;
            durEvents(idx) = [];
            
            if ~isempty(durEvents)
                durEvents(length(durEvents)) = [];
            end
            for i=1:params.binCount
                switch (params.binMethod)
                    case 'log2'
                        bin(o,i) = sum((durEvents >= 2^i) & (durEvents <= 2^(i+1)-1));
                        if (i == params.binCount) && params.includeHighValues
                            bin(o,i) = sum(durEvents >= 2^i);
                        end
                    otherwise
                        bin(o,i) = sum(sum(durEvents == (i-1)*params.binWidth+1:i*params.binWidth));
                        if (i == params.binCount) && params.includeHighValues
                            bin(o,i) = sum(durEvents >= (i-1)*params.binWidth+1);
                        end
                end
            end
        end
        meanVals(:,a) = nanmean(bin);
        dev = nanstd(bin);
        SEM(:,a) = dev/sqrt(numel(bin));
        legendstr(a) = objlist{a}(1).getInfo('Starvation');
    end
    
    if ~isempty(params.saveDataAs)
        %Saves data along with exp ID labels
        cBinVals = arrayfun(@(x) x, bin, 'UniformOutput', false);
        xlsdata = cell(numel(obj), params.binCount+1);
        xlsdata(:,1) = obj.getInfo('SectorID');
        xlsdata(:,2:end) = cBinVals;
        xlswrite(params.saveDataAs,xlsdata);
    end
    
    bar(meanVals);

    ngroups = params.binCount;
    nbars = numel(objlist);
    groupwidth = min(0.8, nbars/(nbars + 1.5));

    for a = 1:numel(objlist)
        xErr = (1:ngroups) - groupwidth/2 + (2*a-1) * groupwidth / (2*nbars);
        errorbar(xErr,meanVals(:,a),SEM(:,a),SEM(:,a),'k','LineStyle','none');
    end

    switch (params.binMethod)
        case 'log2'
            xlabel('Stop duration(log2)','FontSize',params.fontSize,'FontWeight','bold');
        otherwise
            xlabel(['Stop duration(',num2str(params.binWidth),'s/bin)'],...
                'FontSize',params.fontSize,'FontWeight','bold');
    end
    ylabel('Avg number of stop events','FontSize',params.fontSize,'FontWeight','bold')
    legend(legendstr);
    return
end

% -------------------------------------------------------------------------

%MovementDistribution
if strcmp(params.defaultPlot,'MovementDistribution')
    v = nan(params.binCount,length(objlist));
    legendstr = cell(length(objlist),1);
    hold on
    
    for c=1:length(objlist)
        co = objlist{c};
        dc = cell(numel(co),1);
        for o=1:numel(co)
            %if time indices, calculate here
            if ~isempty(params.timeWindow)
                tws = co(o).timeToFrame(params.timeWindow(1));
                twe = co(o).timeToFrame(params.timeWindow(2));
                params.indices = tws:twe;
            end
            %d = co(o).distance(params.indices);
            d = co(o).(params.toPlot)(params.indices);
            d = d / co(o).PIX2MM;
            if ~isempty(params.filterFunction)
                l = co(o).(params.filterFunction)(params.filterParams{1},params.indices);
                dc{o} = d(l);
            else
                dc{o} = d;
            end
        end
        dist = cat(1,dc{:});
        dist(dist < params.keepMinimum) = [];
        dist(isnan(dist)) = [];
        histogram(dist,'BinLimits',params.binLimits,'BinWidth',params.binWidth,...
            'Normalization',params.normalization);
        legendstr(c) = objlist{c}(1).getInfo(params.legendLabel);        
    end

    title('Displacement','Interpreter','none')
    xlabel('Displacement (mm)','FontSize',params.fontSize,'FontWeight','bold');
    ylabel('Probability/interval','FontSize',params.fontSize,'FontWeight','bold');    
    legend(legendstr)
    return
end
%--------------------------------------------------------------------------

% Basic plots
switch (params.toPlot)
    case {'isActive','isActiveStacked','isActiveStackedMM','onArea',...
            'farFromArea','nearArea', 'closeToOrOnArea', 'onAreaAll'}
        patchTop = 1;
    case 'distanceToArea'
        patchTop = 15;
    otherwise
        patchTop = 80;
end

binsize = params.binSize;

dst = PanDataHandleList.DAYSTARTTIME;
det = PanDataHandleList.DAYENDTIME;
fps = PanDataHandleList.FPS;
binsPerDay = fps * 3600 * 24/binsize;
binsPerHour = binsPerDay/24;
pstart = det*binsPerHour-dst*binsPerHour;
theDarkSide = dst+24-det;
pend = theDarkSide*binsPerHour+pstart;

maxBinsList = zeros(length(objlist),1);
for c=1:length(objlist)
    obj = objlist{c};
    maxBinsList(c) = max(obj.numBins(binsize));
    %Skims through the list and identifies the maximum bin count
end
maxBinValue = max(maxBinsList);

for i=1:ceil(maxBinValue/binsPerDay)
    pstart2 = pstart + (i-1)*binsPerDay;
    pend2 = pend + (i-1)*binsPerDay;
    patch([pstart2 pend2 pend2 pstart2],[0 0 patchTop patchTop],[0 0 0 0],[0.7 0.7 0.7]);
end

%Creates a dark patch for night hours
legendstr = cell(length(objlist),1);
offsBinsList = zeros(length(objlist),1);

for c=1:length(objlist)
    obj = objlist{c};
    
    startTimes = obj.startTime;
    maxStart = max(startTimes);
    %Retrieves the latest start time
    offsTimeMax = maxStart - dst;
    offsFramesMax = fps * 3600 * offsTimeMax;
    offsBinsMax = round(offsFramesMax/binsize);
    binsPerDay = fps * 3600 * 24/binsize;
    %maxBins should be calculated from data
    
    maxBins = ceil(offsBinsMax + max(maxBinValue, binsPerDay));
    %Readjusts and plots bins acc to the right start time
    binVals = nan(maxBins,numel(obj));
    
    %Loop over AT    
    if params.multiRounds
        if strcmp(params.toPlot,'isActiveStacked')
            AT = obj(1).ACTIVITYTHRESHOLDS;
        elseif strcmp(params.toPlot,'isActiveStackedMM')
            AT = obj(1).ACTIVITYTHRESHOLDSMM;
        end
        rounds = numel(AT)-1;
    else
        rounds = 1;
    end
    
    %TEST
%     meanVals = nan(rounds,size(binVals,2));
%     stdval = nan(rounds,size(binVals,2));
%     semval = nan(rounds,size(binVals,2));
    meanVals = nan(rounds,size(binVals,1));
    stdval = nan(rounds,size(binVals,1));
    semval = nan(rounds,size(binVals,1));    

    saveData = cell(rounds,1);
    for a=1:rounds
        for o=1:numel(obj)
            offsTime = startTimes(o) - dst;
            offsFrames = fps * 3600 * offsTime;
            offsBins = round(offsFrames/binsize);
            
            if isempty(params.plotParams)
                d = obj(o).(params.toPlot);
            else
                if numel(params.plotParams) == 1
                    d = obj(o).(params.toPlot)(params.plotParams{:});
                else
                    % Works only for 2 params
                    d = obj(o).(params.toPlot)(eval(params.plotParams{1}),...
                        eval(params.plotParams{2}));
                end
            end
            switch (params.binMethod)
                case 'sum'
                    bvTemp = nansum(obj(o).binMe(d,binsize),2);
                case 'mean'
                    bvTemp = nanmean(obj(o).binMe(d,binsize),2);
                otherwise
                    warning('I don''t understand what you mean.');
            end
            % Assign values, array will grow if demanded
            % calculated values are pasted onto empty NaN lists
            binVals(offsBins:offsBins+numel(bvTemp)-1,o) = bvTemp;
        end
        meanValues = nanmean(binVals,2);
        %patch as NaN
        if ~isempty(params.subPlotRange)
            t = nan(numel(meanValues),1);
            %fix max
            nVals = numel(meanValues);
            if max(params.subPlotRange) > nVals
                disp('SubPlotRange outside limits, adjusting to ',num2str(nVals))
                params.subPlotRange(params.subPlotRange > nVals) = nVals;
                params.subPlotRange = unique(params.subPlotRange);
            end
            t(params.subPlotRange) = meanValues(params.subPlotRange);
            meanValues = t;
        end
        
        xErr = numel(meanValues);
        
        dev = nanstd(binVals,0,2);
        SEM = dev/sqrt(size(binVals,2));
        if ~isempty(params.subPlotRange)
            t = nan(numel(SEM),1);
            t(params.subPlotRange) = SEM(params.subPlotRange);
            SEM = t;
        end
        
        saveData{a} = binVals;
        
        %From plotDist...
        if strcmp(params.plotType,'stackedBars')
            meanVals(a,:) = nanmean(binVals,2);
            stdval(a,:) = nanstd(binVals,0,2);
            semval(a,:) = stdval(a,:)/sqrt(size(binVals,2));
            xErr = size(meanVals,2);
        end
        % case between
    end
    
    sprintf('Size meanValues %s',size(meanValues));
    sprintf('Size SEM %s',size(meanValues));
    
    switch (params.plotType)
        case 'shadedErrorBar'
            hs = shadedErrorBar(1:xErr,meanValues,SEM,'lineProps',...
                {'Color',params.colorList{c},'LineWidth',4});
            g(c) = hs.mainLine;
        case 'bar'
            g(c) = bar(meanValues);
            errorbar(1:xErr,meanValues,SEM,SEM);
        case 'stackedBars'
            if ~isempty(params.subPlotRange)
                t = nan(size(meanVals));
                t(:,params.subPlotRange) = meanVals(:,params.subPlotRange);
                meanVals = t;
                t = nan(size(semval));
                t(:,params.subPlotRange) = semval(:,params.subPlotRange);
                semval = t;
            end
            g = bar(meanVals','Stacked');
            yPos = cumsum(meanVals,'omitnan');
            for a = 1:numel(AT)-1
                errorbar(1:xErr,yPos(a,:),semval(a,:),semval(a,:),'k','LineStyle','none');
            end
        otherwise
            disp('not implemented yet')
    end
    
    %if ~isempty(params.legendLabel)
    if isfield(params,'legendLabel')
        legendstr{c} = [cell2mat(objlist{c}(1).getInfo(params.legendLabel)),' [',num2str(numel(objlist{c})),']'];
    end    
    offsBinsList(c) = offsBinsMax;
    
    %Save the data
    if ~isempty(params.saveDataAs)
        switch (params.toPlot)
            case {'isActiveStacked','isActiveStackedMM'}
                nAT = numel(AT)-1;
                binVals = [saveData{:}]';
                size(binVals)
                cBinVals = arrayfun(@(x) x, binVals, 'UniformOutput', false);
                xlsdata = cell(numel(obj)*nAT, size(binVals,2)+1);                
                expIDs = obj.getInfo('ExptID');
                for a=1:nAT
                    ids(numel(obj)*(a-1)+1:numel(obj)*a) = strcat(expIDs,'_',num2str(AT(a)),'_',num2str(AT(a+1)));
                end
                xlsdata(:,1) = ids;
                xlsdata(:,2:end) = cBinVals;
            otherwise
                cBinVals = arrayfun(@(x) x, binVals, 'UniformOutput', false);
                xlsdata = cell(numel(obj), maxBins+1);
                xlsdata(:,1) = obj.getInfo('SectorID');
                xlsdata(:,2:end) = cBinVals';
        end
        xlswrite(params.saveDataAs,xlsdata);
    end
end

titlestr = strcat(objlist{1}(1).getInfo('Genotype'));

switch (params.xTicks)
    %default number of ticks plotted if otherwise indicated
    case 'default'
        maxTicks = (maxBinValue/2/binsPerHour)*(fps * 3600 * 2/binsize)+max(offsBinsList);
        xticks(fps * 3600/binsize:fps * 3600 * 2/binsize:maxTicks)
        xl = dst+1:2:(maxTicks+dst+1)*2;
    case '12h'
        maxTicks = (maxBinValue/6/binsPerHour)*(fps * 3600 * 6/binsize)+max(offsBinsList);
        xticks(fps * 3600/binsize:fps * 3600 * 6/binsize:maxTicks)
        xl = dst+1:12:(maxTicks+dst+1)*2;
    otherwise
        if ~isnumeric(params.xTicks)
            maxTicks = (maxBinValue/2/binsPerHour)*(fps * 3600 * 2/binsize)+max(offsBinsList);
            xticks(fps * 3600/binsize:fps * 3600 * 2/binsize:maxTicks)
            xl = dst+1:2:(maxTicks+dst+1)*2;
        else
            %To plot any chosen number of ticks by entering the numeric value
            maxTicks = (maxBinValue/params.xTicks/binsPerHour)*(fps * 3600 * params.xTicks/binsize)+max(offsBinsList);
            xticks(fps * 3600/binsize:fps * 3600 * params.xTicks/binsize:maxTicks)
            xl = dst+1:params.xTicks:(maxTicks+dst+1)*2;
        end
end

for i=1:ceil(maxBinValue/binsPerDay)
    idx = xl > 24;
    xl(idx) = xl(idx)-24;
end
xticklabels(xl);

h.Children.XAxis.FontSize = params.fontSize;
h.Children.YAxis.FontSize = params.fontSize;

xlabel('Time of Day (hr)','FontSize',params.fontSize,'FontWeight','bold');

switch (params.toPlot)
    case 'isActive'
        h.Name = ['Walking Activity for ',titlestr{1}];
        title(['Walking Activity for ', titlestr],'Interpreter','none')
        ylabel('Fraction of active flies','FontSize',params.fontSize,'FontWeight','bold');
    case 'onArea'
        h.Name = ['Fraction of Time on Patch for ',titlestr{1}];
        title(['Fraction of Time on Patch for ', titlestr],'Interpreter','none')
        ylabel('Fraction of time spent on food patch','FontSize',params.fontSize,'FontWeight','bold');
    case 'onAreaAll'
        h.Name = ['Fraction of Time on Patch for ',titlestr{1}];
        title(['Fraction of Time on Patch for ', titlestr],'Interpreter','none')
        ylabel('Fraction of time spent on food patch','FontSize',params.fontSize,'FontWeight','bold');
    case 'farFromArea'
        h.Name = ['Fraction of Time away from the Patch for ',titlestr{1}];
        title(['Fraction of Time away from the Patch for ', titlestr],'Interpreter','none')
        ylabel('Fraction of time spent on food patch','FontSize',params.fontSize,'FontWeight','bold');
    case 'nearArea'
        h.Name = ['Fraction of Time near the Patch for ',titlestr{1}];
        title(['Fraction of Time near the Patch for ', titlestr],'Interpreter','none')
        ylabel('Fraction of time spent on food patch','FontSize',params.fontSize,'FontWeight','bold');
    case 'closeToOrOnArea'
        h.Name = ['Fraction of Time near the Patch for ',titlestr{1}];
        title(['Fraction of Time near the Patch for ', titlestr],'Interpreter','none')
        ylabel('Fraction of time spent on food patch','FontSize',params.fontSize,'FontWeight','bold');
    case 'distanceToArea'
        h.Name = ['Distance from Patch for ',titlestr{1}];
        title(['Distance from Patch for ', titlestr],'Interpreter','none')
        ylabel('Average distance (mm)','FontSize',params.fontSize,'FontWeight','bold');
    case 'isActiveStacked'
        h.Name = ['Change ',titlestr{1}];
        title(['Walking Activity for ', titlestr],'Interpreter','none')
        ylabel('Fraction of active flies','FontSize',params.fontSize,'FontWeight','bold');
        for i=1:a
            legendstr(i) = {['Text [',num2str(AT(i)),'-',num2str(AT(i+1)),']']};
        end
    case 'isActiveStackedMM'
        h.Name = ['Change ',titlestr{1}];
        title(['Walking Activity for ', titlestr],'Interpreter','none')
        ylabel('Fraction of active flies','FontSize',params.fontSize,'FontWeight','bold');
        for i=1:a
            legendstr(i) = {[num2str(AT(i)),' - ',num2str(AT(i+1)),' mm']};
        end
    otherwise
        h.Name = ['Displacement for ',titlestr{1}];
        title(['Displacement for ', titlestr],'Interpreter','none')
        ylabel('Average distance (mm)','FontSize',params.fontSize,'FontWeight','bold');
end

axis ([0 maxTicks 0 patchTop]);
legend(g, legendstr)
end
