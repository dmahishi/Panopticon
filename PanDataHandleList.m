classdef PanDataHandleList < BasicDataHandleList
%PANDATAHANDLELIST Class to handle Panopticon experiment data
%   PanDataHandleList is a subclass of BasicDataHandleList specifically
%   designed to handle data generated in the Panopticon experiment
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
    
    properties
        areas % areas (arena and foodpatch if available)
        border % points defining sectors
    end
    
    properties (Constant = true, Hidden = true)
        NVALS = 1800; % default expected batch size
        IMGRES = 1024; % image resolution
        FPS = 1; % framerate
        SECTCOUNT = 8; % number of sectors
        DAYSTARTTIME = 7; % lights on time (subjective)
        DAYENDTIME = 21; % lights off time (subjective)
        ACTIVITYTHRESHOLDS = [2 30 100 300]; % Movement thresholds in Pixel
        ACTIVITYTHRESHOLDSMM = [0 2 13 30]; %[0.2 2 30], in Millimeter
        RESSCALE = 1/8; % Scaling fr residency plot
        IMGFILETYPE = '.jpeg';
        PIX2MM = 978/83; % about 12 pix/mm
    end
    
    properties (Hidden = true)
        defVal % default value (i.e. all NaNs)
    end
    
    methods
        %constructor
        function obj = PanDataHandleList(data)
            if nargin > 0
                obj.data = data;
            end
            obj.defVal.Area = nan;
            obj.defVal.x = nan;
            obj.defVal.y = nan;
            obj.defVal.Slice = nan;
            obj.defVal.Sector = nan;
            obj.defVal.Circ = nan;
        end
        
        function dhl = sLIQ(obj,varargin)
            %To sort data according to specific field names from metadata
            dhl = obj.selectListsByInfoQuery(varargin{:});
        end
        
        %Short utility functions
        function shiftCenter(obj,x,y)
            %to shift the center of the arena
            for o=1:numel(obj)
                obj(o).areas(1).center.x = obj(o).areas(1).center.x + x;
                obj(o).areas(1).center.y = obj(o).areas(1).center.y + y;
            end
        end

        function setAreas(obj,a)
            for o=1:numel(obj)
                obj(o).areas = a;
            end
        end
        
        function addArea(obj,paramin)
            if nargin < 2
                paramin = {};
            end
            p = inputParser;
            p.KeepUnmatched = 1;
            p.parse(paramin{:});
            params = p.Unmatched;            
            
            if isfield(params,'id')
               % to add id for the new area(s)
                id = params.id;
            else
                id = numel(obj.areas)+1;
            end
            
            if isfield(params,'x')
               % to add x position (cartesian)
                obj.areas(id).center.x = params.x;
            end
            
            if isfield(params,'y')
               % to add y position (cartesian)
                obj.areas(id).center.y = params.y;
            end
                        
            if isfield(params,'radius')
                % to add y position (cartesian)
                obj.areas(id).radius = params.radius;
            end
            
            if isfield(params,'name')
                obj.areas(id).name = params.name;
            else
                obj.areas(id).name = 'VirtualPatch';
            end
            
        end
        
        function makeVirtualPatches(obj,paramin)
            % To draw virtual food spots or patches on the Arena
            if nargin < 2
                paramin = {};
            end
            p = inputParser;
            p.KeepUnmatched = 1;
            p.parse(paramin{:});
            params = p.Unmatched;
            
            if ~isfield(params,'id')
               % id to the virtual patch(es) (2 or greater, 1 = Arena)
                params.id = numel(obj(1).areas)+1;
            end
            
            if ~isfield(params,'offsetAngle')
               % reference angle for virtual patch(es) position(s)             
                params.offsetAngle = 0.2;
            end
            
            if ~isfield(params,'offsetDist')
               % reference dist for virtual patch(es) position(s)
                params.offsetDist = 1.35;
            end
            
            if ~isfield(params,'radius')
               % radius measure for virtual patch(es) position(s)
                params.radius = 20;
            end
            
            if ~isfield(params,'name')
               % reference name(s) for virtual patch(es) position(s)
                params.name = 'VirtualPatch';
            end
            
            for o=1:numel(obj)
                sector = cell2mat(obj(o).getInfo('Sector'));
                x = mean([obj(o).border(sector,1) ...
                    obj(o).border(sector+1,1) ...
                    obj(o).areas(1).center.x]);
                y = mean([obj(o).border(sector,2) ...
                    obj(o).border(sector+1,2) ...
                    obj(o).areas(1).center.y]);
                xloc = x - obj(o).areas(1).center.x;
                yloc = y - obj(o).areas(1).center.y;
                [theta, rho] = cart2pol(xloc,yloc);
                if params.id == 2
                    [xshift, yshift] = pol2cart(theta,rho*1.1);
                elseif params.id >= 3
                    [xshift, yshift] = pol2cart(theta+params.offsetAngle,rho*params.offsetDist);
                end
                xshift = xshift + obj(o).areas(1).center.x;
                yshift = yshift + obj(o).areas(1).center.y;
                obj(o).addArea({'x', xshift, 'y', yshift, 'radius', params.radius, ...
                    'id', params.id, 'name', params.name})
            end
        end
            
        function [theta, rho] = toPolar(obj)
            % to measure and set polar co-ordinates
            xloc = obj.get('x') - obj.areas(1).center.x;
            yloc = obj.get('y') - obj.areas(1).center.y;
            [theta, rho] = cart2pol(xloc,yloc);
        end
        
        function d = distanceToArea(obj,nr,indices)
            %Calculates distance from specific co-ordinates
            if nargin < 3
                indices = ':';
            end
            xloc = obj.get('x',indices) - obj.areas(nr).center.x;
            yloc = obj.get('y',indices) - obj.areas(nr).center.y;
            [~, rho] = cart2pol(xloc,yloc);
            d = rho;
        end
        
        function l = onArea(obj,nr,indices)
            %fraction of time spent on the food patch
            if nargin < 3
                indices = ':';
            end
            xloc = obj.get('x',indices) - obj.areas(nr).center.x;
            yloc = obj.get('y',indices) - obj.areas(nr).center.y;
            [~, rho] = cart2pol(xloc,yloc);
            l = rho <= obj.areas(nr).radius;
        end
        
        function l = onAreaAll(obj,nr,indices)
            %fraction of time spent on the food patch
            if nargin < 3
                indices = ':';
            end
            
            l = false(sum(obj.getLength),1);            
            idx = [0; cumsum(obj.getLength)];
            
            for o=1:numel(obj)
                x = obj(o).get('x',indices) - obj(o).areas(nr).center.x;
                y = obj(o).get('y',indices) - obj(o).areas(nr).center.y;
                [~, rho] = cart2pol(x,y);
                l(idx(o)+1:idx(o+1)) = rho <= obj(o).areas(nr).radius;
            end            
        end        
        
        function l = closeToOrOnArea(obj,nr,radius,indices)
            %fraction of time around and on the food patch
            if nargin < 4
                indices = ':';
            end
            xloc = obj.get('x',indices) - obj.areas(nr).center.x;
            yloc = obj.get('y',indices) - obj.areas(nr).center.y;
            [~, rho] = cart2pol(xloc,yloc);
            l = rho <= obj.areas(nr).radius+radius;
        end
        
        function l = farFromArea(obj,nr,indices)
            %fraction of time away from the food patch
            if nargin < 3
                indices = ':';
            end
            xloc = obj.get('x',indices) - obj.areas(nr).center.x;
            yloc = obj.get('y',indices) - obj.areas(nr).center.y;
            [~, rho] = cart2pol(xloc,yloc);
            l = rho > obj.areas(nr).radius+40;
        end
        
        function l = nearArea(obj,nr,indices)
            %fraction of time around the food patch but not on
            if nargin < 3
                indices = ':';
            end
            xloc = obj.get('x',indices) - obj.areas(nr).center.x;
            yloc = obj.get('y',indices) - obj.areas(nr).center.y;
            [~, rho] = cart2pol(xloc,yloc);
            l = (rho > obj.areas(nr).radius) & (rho < obj.areas(nr).radius+20);
        end
        
        function f = fractionOnArea(obj,nr,indices)
            %Percentage of time flies spend on the patches
            if nargin < 3
                indices = ':';
            end
            f = nan(numel(obj),1);
            for o=1:numel(obj)
                f(o) = nanmean(obj(o).onArea(nr,indices));
            end
        end
        
        function d = distance(obj,indices)
            %for calculation of diaplacement in pixel units
            if nargin < 2
                indices = ':';
            end
            x = obj.get('x',indices);
            y = obj.get('y',indices);

            xd = diff(x);
            yd = diff(y);
            d = hypot(xd,yd);
            d(end+1) = nan;
        end
        
        function d = distanceInMM(obj,indices)
            %for calculation of diaplacement in mm units
            if nargin < 2
                indices = ':';
            end
            d = obj.distance(indices) / obj.PIX2MM;
        end
        
        
        function d = isActiveBetween(obj,minThres,maxThres,indices)
            %Percentage active above and below each threshold excluding NaNs
            if nargin < 4
                indices = ':';
            end
            %activity status of each frame
            dist = obj.distance(indices);
            act = (dist > minThres) & (dist <= maxThres);
            act_d = double(act);
            act_d(isnan(obj.distance)) = nan;
            d = act_d;
        end

        function d = isActiveBetweenMM(obj,minThres,maxThres,indices)
            %Percentage active above and below each threshold excluding NaNs
            if nargin < 4
                indices = ':';
            end
            %activity status of each frame
            dist = obj.distance(indices);
            dist = dist / obj.PIX2MM;
            act = (dist > minThres) & (dist <= maxThres);
            act_d = double(act);
            act_d(isnan(obj.distance)) = nan;
            d = act_d;
        end
        
        function d = isActive(obj,indices)
            % Frames with activity above threshold
            if nargin < 2
                indices = ':';
            end
            minThres = obj.ACTIVITYTHRESHOLDS(1);
            %activity status of each frame
            dist = obj.distance(indices) > minThres;
            dist_d = double(dist);
            dist_d(isnan(obj.distance(indices))) = nan;
            d = dist_d;
        end
        
        function f = fractionActive(obj,indices)
            %Calculates the percentage of flies active
            if nargin < 2
                indices = ':';
            end
            f = nan(numel(obj),1);
            for o=1:numel(obj)
                f(o) = nanmean(obj(o).isActive(indices));
            end
        end
        
        function [posEventStarts, durEvents] = calculateStops(obj,indices)
            %Calculates stop events from across the entire data set
            if nargin < 2
                indices = ':';
            end

            act = obj.isActive(indices);
            act(isnan(act)) = 4;
            da = diff(act);
            da = [nan; da];
            %Find the beginning of stop events
            posEventStarts = find((da==-1) + (da==-4));
            %Find the end of stop events
            posEventEnds = find((da==1) + (da==4));
            if act(1) == 0
                posEventEnds(1) = [];
            end
            if (numel(posEventStarts) > numel(posEventEnds))
                posEventEnds(numel(posEventStarts)) = numel(act);
            end
            durEvents = posEventEnds - posEventStarts;
        end
        
        function bd = binMe(~,input,binsize)
            %For creating bins
            dl = ceil(numel(input)/binsize);
            bd = nan(dl,binsize);
            for i=1:dl
                %Start of bin
                si = (i-1)*binsize+1;
                %End of bin
                ei = i*binsize;
                ta = nan(binsize,1);
                %Creates a nan list of binsize length
                ti = input(si:min(ei,numel(input)));
                %Takes upto all remaining values
                ta(1:length(ti)) = ti;
                %Puts the above values into ta
                bd(i,:) = ta;
                %Assigns to bin(s)
            end
        end
        
        function d = numBins(obj,binsize)
            %total number of bins
            d = (obj.getLength / binsize);
        end
        
        function d = maxBins(obj,binsize)
            %Maximum number of bins
            d = max(obj.numBins(binsize));
        end
        
        function assignSector(obj)            
            % Assigns Sectors for fly positions using polar coordinates            
            for o=1:numel(obj)
                xloc = obj(o).border(:,1) - obj(o).areas(1).center.x;
                yloc = obj(o).border(:,2) - obj(o).areas(1).center.y;
                [theta, rho] = cart2pol(xloc,yloc);
                [xshift, yshift] = pol2cart(theta,rho*1.1);
                xshift = xshift + obj(o).areas(1).center.x;
                yshift = yshift + obj(o).areas(1).center.y;
                x = obj(o).get('x');
                y = obj(o).get('y');
                sectors = nan(numel(x),1);
                for i=1:size(obj(o).border,1)-1
                    insect = inpolygon(x,y,[xshift(i),xshift(i+1),obj(o).areas(1).center.x],...
                        [yshift(i),yshift(i+1),obj(o).areas(1).center.y]);
                    sectors(insect) = i;
                end
                obj(o).addFieldAndValues('Sector',sectors);
            end
        end
        
        function addImagePath(obj)
            %adds the path of the image-containing folder
            for o=1:numel(obj)
                dataPath = obj(o).getInfo('DataPath');
                files = dir([dataPath{:},'\EXP_*']);
                if files(1).isdir
                    obj(o).addInfo('ImagePath',cell2mat(strcat(dataPath,'\',files(1).name)));
                else
                    warning(['No Image folder found under ',dataPath]);
                end
            end
        end
        
        function addDataPath(obj, datapath)
            %adds the path of the image-containing folder
            for o=1:numel(obj)
                if (exist(datapath, 'dir') == 7)
                    obj(o).addInfo('DataPath',datapath);
                else
                    warning(['Folder not found ',datapath]);
                end
            end
        end
        
        function s = findImageForFrame(obj,idx)
            %Finds the corresponding image by frame number idx
            path = cell2mat(obj(1).getInfo('ImagePath'));
            files = dir([path '\*' obj(1).IMGFILETYPE]);
            s = [path '\' files(idx).name];
        end        
        
        function h = showLongDistanceMoves(obj,nFrames)
            %Finding frames which have the longest displacement measured
            if nargin < 2
                nFrames = 10;
            end
            for o=1:numel(obj)
                [~,idx] = maxk(obj(o).distance,nFrames);
                for i=1:numel(idx)
                    h(i) = obj.plotTracksOnImage(idx(i)-1:idx(i)+1); %#ok<AGROW>
                end
            end
        end
        
        function h = plotTracksOnImage(obj,idx)
            % load and show image corresponding to frame
            % Sanity check function, only a single folder will be looked at
            path = cell2mat(obj(1).getInfo('ImagePath'));
            files = dir([path '\*' obj(1).IMGFILETYPE]);
            h = figure;
            for i=1:numel(idx)
                imgs(:,:,:,i) = imread([path '\' files(idx(i)).name]);
            end
            img = min(imgs,[],4);
            imshow(img)
            hold on
            for o=1:numel(obj)
                plot(obj(o).get('x',idx),obj(o).get('y',idx),'o-','MarkerSize',10);
            end
            h.Name = [files(idx(1)).name,' is frame ',num2str(idx(1))];
            set(gca, 'ColorOrder', parula(obj(1).SECTCOUNT));
        end
        
        function h = plotXY(obj,indices)
            %To plot the position of food patches, for specific time bins
            if nargin < 2
                indices = ':';
            end
            size = 20;
            h = figure;
            hold on
            for o=1:numel(obj)
                scatter(obj(o).get('x',indices),obj(o).get('y',indices),size);
                c = [obj(o).areas.center];
                x = [c.x]';
                y = [c.y]';
                r = [obj(o).areas.radius];
                viscircles([x y],r);
            end
            axis([0 obj(o).IMGRES 0 obj(o).IMGRES])
            axis square
            set(gca,'Ydir','reverse')
            set(gca, 'ColorOrder', parula(obj(1).SECTCOUNT));
        end
        
        function h = plotWalkingTrace(obj,indices)
            %to plot the walking trail of the fly as a fun of time
            if nargin < 2
                indices = ':';
            end
            h = figure;
            hold on
            for o=1:numel(obj)
                x = obj(o).get('x',indices);
                y = obj(o).get('y',indices);
                c = 1:numel(x);
                surface([x,x], [y,y], [c(:),c(:)], 'EdgeColor','flat', 'FaceColor','none');
                colormap(jet(numel(x)));
                set(gca,'Ydir','reverse');
                zlabel('Time');
            end
        end
        
        function h = plotPositions(obj,indices)
            %to plot the position of the fly as a fun of time
            if nargin < 2
                indices = ':';
            end
            h = figure;
            hold on
            for o=1:numel(obj)
                x = obj(o).get('x',indices);
                y = obj(o).get('y',indices);
                c = 1:numel(x);
                scatter3(x(:), y(:), c(:),10,c(:),'filled');
                colormap(jet(numel(x)));
                set(gca,'Ydir','reverse');
                zlabel('Time');
            end
        end
        
        function h = plotResidency(obj)
            %To plot the location probability of flies as a heatmap
            if nargin == 1
                t = obj.flattenLists;
            end
            h = figure;
            hold on
            x = t.get('x');
            y = t.get('y');
            res = t.IMGRES*t.RESSCALE;
            image(hist3([[0;x(:);t.IMGRES], [0;y(:);t.IMGRES]], [res res])');
            colormap(hot);
            h.Children.YDir = 'reverse';
            axis square
        end
        
        function l = isNan(obj,indices)
            %Returns NaNs
            if nargin < 2
                indices = ':';
            end
            l = isnan(obj.get('x',indices));
        end
        
        function d = fractionNan(obj,indices)
            %To calculate the percentage of NaNs in a dataset
            if nargin < 2
                indices = ':';
            end
            d = nan(numel(obj),1);
            for o=1:numel(obj)
                d(o) = mean(obj(o).isNan(indices));
            end
        end
              
        function d = startTime(obj)
            %The actual start time of all expts as entered in metaData
            d = nan(numel(obj),1);
            for o = 1:numel(obj)
                strTime = num2str(cell2mat(obj(o).getInfo('StartTime')),'%06d');
                d(o) = str2num(strTime(1:2))+str2num(strTime(3:4))/60 + str2num(strTime(5:6))/3600;
            end
        end
        
        function d = endTime(obj)
            %The endtime calculated according to the length of the data
            d = nan(numel(obj),1);
            for o = 1:numel(obj)
                d(o) = obj(o).frameToTime(obj(o).getLength);
            end
        end
        
        function d = frameToTime(obj,frames,to24h)
            %Coversion of Frames to actual time of the day
            if nargin < 3
                to24h = false;
            end
            d = nan(numel(obj),1);
            for o=1:numel(obj)
                t = (frames / (obj(o).FPS * 3600)) + obj(o).startTime;
                if to24h
                    t(t > 24) = t(t > 24) - 24;
                end
                d(o) = t;
            end
        end
        
        function d = timeToFrame(obj,times)
            %Conversion of time of the day to frame numbers
            d = nan(numel(obj),1);
            for o=1:numel(obj)
                d(o) = int32((times - obj(o).startTime) * obj(o).FPS * 3600);
            end
        end
        
        function pdhl = syncStartTimes(obj)
            % To synchronise the start time to the earliest hr in a list
            % for plotting. Returns a new list with shifted start times
            pdhl = obj.copy;
            st = pdhl.getInfo('StartTime');
            min_st = min([st{:}]);
            for o=1:numel(pdhl)
                pdhl(o).addInfo('StartTime',min_st);
            end
        end
        
        function h = plotNans(obj, paramin)
            %to plot the fraction of NaNs (i.e. missed detections)
            if nargin < 2
                paramin = {};
            end
            p = inputParser;
            p.KeepUnmatched = 1;
            p.parse(paramin{:});
            params = p.Unmatched;
            
            if isfield(params,'indices')
                indices = params.indices;
            else
                indices = ':';
            end
            
            h = figure;
            num = numel(obj);
            for o=1:num
                subplot('position',[0.05, 1-o/(num+2)-0.03, 0.75, 1/(num+10)]);
                image(isnan(obj(o).get('x', indices))');
                yticklabels('');
                colormap([1 1 1; 0 0 0])
                fn = max(obj(o).fractionNan(indices),10^-10);
                subplot('position',[0.85, 1-o/(num+2)-0.03, 0.2, 1/(num+10)]);
                pie(fn);
                colormap([1 1 1; 0.4 0.75 0.5])
            end
        end
        
        function h = plotOnArea(obj,paramin)
            %to plot time distribution and fraction on area
            if nargin < 2
                paramin = {};
            end
            p = inputParser;
            p.KeepUnmatched = 1;
            p.parse(paramin{:});
            params = p.Unmatched;
            
            if isfield(params,'indices')
                indices = params.indices;
            else
                indices = ':';
            end
            
            if isfield(params,'target')
                target = params.target;
            else
                target = 2;
            end
            
            if isfield(params,'sorted')
                sorted = params.sorted;
            else
                sorted = false;
            end
            
            h = figure;
            num = numel(obj);
            if sorted
                idx = nan(numel(obj),1);
                for o=1:numel(obj)
                    f = find(obj(o).onArea(target,1:obj(o).NVALS),1,'first');
                    if ~isempty(f)
                        idx(o) = f;
                    end
                end
                [~,sidx] = sort(idx);
            else
                sidx = 1:num;
            end
            for o=1:num
                subplot('position',[0.05, 1-o/(num+2)-0.03, 0.75, 1/(num+10)])
                image(obj(sidx(o)).onArea(target, indices)')
                yticklabels('');
                colormap([1 1 1; 0 0 0])
                ft = max(obj(sidx(o)).fractionOnArea(target, indices),10^-10);
                subplot('position',[0.85, 1-o/(num+2)-0.03, 0.2, 1/(num+10)])
                pie(ft)
                colormap([1 1 1; 0.4 0.75 0.5])
            end            
        end
        
        function repairData(obj,targetSize)
            % Handle duplicate and missed detections
            for o=1:numel(obj)
                ymat = 1:targetSize;
                xmat = obj(o).get('Slice');
                area = obj(o).get('Area');
                matrix = false(length(xmat),length(ymat));
                for i = 1:length(ymat)
                    matrix(:,i) = ismember(xmat,ymat(i));
                end
                sv = sum(matrix);
                badval = find(sv > 1);
                delval = [];
                for i = 1:numel(badval)
                    idx = xmat==badval(i);
                    areaval = area(idx);
                    fidx = find(idx);
                    maxval = max(areaval);
                    delval = [delval; fidx((areaval < maxval))];
                end
                obj(o).data(delval) = [];
                for i=targetSize:-1:1
                    td(i) = obj(o).defVal;
                    td(i).Slice = i;
                    td(i).Sector = o;
                end
                snr = obj(o).get('Slice');
                td(snr) = obj(o).data;
                obj(o).data = td;
            end
        end
        
        function dhls = splitAndRepair(obj,fieldname,targetSize)
            % repairs our data by sorting the right values for each frame/slice(works for Sectors only)
            % Creates an empty list for the mentioned fieldname
            allSect = [];
            %Check if last file was complete (i.e. >= 99%)
            if max(obj(end).get('Slice')) < 0.99*obj(end).NVALS
                obj(end) = [];
            end
            for o=1:numel(obj)
                % Adds values to the above empty list
                allSect = [allSect; obj(o).get(fieldname)];
            end
            uf = unique(allSect);
            for o=1:numel(obj)
                disp(['Splitting and repairing list ' num2str(o)])
                dhl = obj(o).splitByFieldName(fieldname);
                % did not create a list for each expected sector
                if numel(uf) ~= numel(dhl)
                    tdhl = PanDataHandleList.empty(0,8);
                    for i=1:numel(uf)
                        for j=1:numel(dhl)
                            % compare unique sector to sector found in list
                            if uf(i) == dhl(j).get(fieldname,1)
                                tdhl(i) = dhl(j);
                            end
                        end
                        if numel(tdhl) < i
                            for z=obj(1).NVALS:-1:1
                                nanData(z) = obj(1).defVal;
                                nanData(z).Slice = z;
                                nanData(z).Sector = i;
                            end
                            tdhl(i) = PanDataHandleList(nanData);
                        end
                    end
                    dhl = tdhl;
                end
                dhl.repairData(targetSize);
                dhlsT(o,:) = dhl;
            end
            for i = 1:size(dhlsT,2)
                dhls(i) = dhlsT(:,i).flattenLists;
            end
            switch (obj(1).info.ExperimentType)
                case 'Basic'
                    for i=1:numel(dhls)
                        dhls(i).areas = obj(1).areas;
                        dhls(i).border = obj(1).border;
                    end
                case 'Foraging'
                    for i=1:numel(dhls)
                        dhls(i).areas = obj(1).areas(1);
                        sectors = [obj(1).areas.Sector];
                        dhls(i).areas(2) = obj(1).areas(find(sectors == i));
                        dhls(i).border = obj(1).border;
                    end
            end           
        end
        
        function fixGaps(obj)
            % Repair gaps of one single NaNs between frames
            for o=1:numel(obj)
                ln = obj(o).isNan;
                ln(1) = 0; %error handling
                dn = diff(ln);
                gs = find(dn == 1)+1; %gap start
                ge = find(dn == -1)+1; %gap end
                gl = ge - gs; %gap length
                gp = gs(gl == 1); %check if +/-1 out of border
                x = num2cell(mean([obj(o).get('x',gp-1) obj(o).get('x',gp+1)],2));
                y = num2cell(mean([obj(o).get('y',gp-1) obj(o).get('y',gp+1)],2));
                [obj(o).data(gp).x] = deal(x{:});
                [obj(o).data(gp).y] = deal(y{:});
                obj(o).info.fixGaps = 1;
            end            
        end
        
        function db = addMetaData(obj,filename)
            %Loads the Metadata xlsx file and applies field info
            dhl = obj(1).copy;
            db = dhl.loadMetaData(filename);
            if numel(obj) ~= (size(db.Pano,1)-2)
                error('This won''t work, please adjust xlsx file (nr. of rows)')
            end
            fns = db.Pano(1,:);
            for o=1:numel(obj)
                for f = 1:numel(fns)
                    v = db.Pano(o+2,f);
                    obj(o).addInfo(fns{f},v{1});
                end
            end
        end
        
    end
    
    methods %override
        
        function saveToFile(obj,filepath)
            %Saves data to the designated filepath
            disp('Saving, please wait...')
            data = cell(numel(obj),1);
            for o=1:numel(obj)
                data{o} = obj(o).getData;
            end
            info = obj.getInfo;
            area = {obj.areas};
            save(filepath,'data','info','area');
        end
        
    end
    
    methods (Static)
        function dhl = CreateFromFile(filepath,sep)
            %Prepares data for analysis and defines sectors
            if ~exist(filepath,'file')
                error('Could not find file.')
            end
            if endsWith(filepath,'\')
                error('Please do not include the last character!')
            end
            if nargin < 2
                sep = '\t';
            end
            fprintf('Loading file %s\n',filepath)
            csvdata = dlmread(filepath,sep,1,0);
            
            for i=length(csvdata):-1:1
                data(i).Area = csvdata(i,2);
                data(i).x = csvdata(i,3);
                data(i).y = csvdata(i,4);
                data(i).Slice = csvdata(i,6);
                data(i).Circ = csvdata(i,5);
            end
            dhl = PanDataHandleList(data);
            dhl.info.file = filepath;
        end
        
        function dhls = CreateFromFolder(folderpath,sep)
            %Loading data file from the indicated filepath
            %For assigning arena centre and radii for food Patches too
            if ~exist(folderpath,'dir')
                error('Could not find folder.');
            end
            if nargin < 2
                sep = '\t';
            end
            files = dir([folderpath '\*.xls']);
            csvfilename = [folderpath '\Centre.csv'];
            if ~exist(csvfilename,'file')
                warning('csvfile not found!')
                return
            end
            csvfile = csvread(csvfilename,1,0);
            center.x = csvfile(:,3);
            center.y = csvfile(:,4);
            
            radii = round(sqrt(csvfile(:,2)/pi));
            areas(1).name = 'Arena';
            areas(1).center.x = center.x(1);
            areas(1).center.y = center.y(1);
            areas(1).radius = radii(1);
            areas(1).Sector = 0;

            borderFile = [folderpath '\BorderPoints.csv'];
            if exist(borderFile,'file')
                borderPoints = csvread(borderFile,1,0);
                borderPoints(9,:) = borderPoints(1,:);
            else
                warning('BorderPoints.csv not found!')
            end
            
            %decide if there are food spots
            if numel(radii) == 1
                ExperimentType = 'Basic';
            else
                ExperimentType = 'Foraging';
                for i=numel(radii):-1:2                    
                    areas(i).name = 'Foodpatch';
                    areas(i).center.x = center.x(i);
                    areas(i).center.y = center.y(i);
                    areas(i).radius = radii(i);
                    %Set Sector Info
                    areas(i).Sector = i-1;
                end
            end
            
            for f=1:numel(files)
                dhls(f) = PanDataHandleList.CreateFromFile([files(f).folder,'\',files(f).name],sep);
                dhls(f).setAreas(areas);
                dhls(f).info.ExperimentType = ExperimentType;
                dhls(f).border = borderPoints(:,3:4);
            end
            
        end
        
        function dhls = LoadFromFile(filepath)
            %To load specific .mat(or other) files using the filepath
            disp('Loading, please wait...')
            data = load(filepath,'data');
            info = load(filepath,'info');
            area = load(filepath,'area');
            for o = numel(data.data):-1:1
                dhls(o) = PanDataHandleList(data.data{o});
                dhls(o).info = info.info(o);
                if isfield(area,'area')
                    dhls(o).areas = area.area{o};
                end
            end
        end
        
        function dhls = LoadFromTextFile(filepath)
            %To load specific .mat(or other) files using the filepath
            disp('Loading experiments, please wait...')
            if exist(filepath,'file') == 2
                fid = fopen(filepath,'r');
            else
                error(['File ',filepath,' not found'])
            end
            C = textscan(fid, '%s');
            files = C{:};
            dhl = PanDataHandleList.empty(numel(files),0);
            for i=1:numel(files)
                fn = files(i);
                dhl{i} = PanDataHandleList.LoadFromFile(fn{:});
            end
            dhls = [dhl{:}];
        end
        
    end
    
end
