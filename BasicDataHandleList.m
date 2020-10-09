classdef BasicDataHandleList < matlab.mixin.Copyable
%BASICDATAHANDLELIST Basic class to handle experiment data for structs
%   This class can handle generic data backed by structs and/or arrays of 
%   structs. Specific metadata can be attached to each object. Create 
%   specific subclasses for detailed analysis.
%
%   Author: Tilman Triphan, tilman.triphan@uni-leipzig.de
%   License: GNU General Public License v3.0
%   Copyright (C) 2020  Tilman Triphan
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
    
    properties (SetAccess = protected, GetAccess = protected)
        data; %underlying experiment data
        info; %metadata
        hash; %hash, used for (some) fast(er) computations
        db;   %database (i.e. xlsx table) for metadata
    end
    
    methods        
        
        function obj = BasicDataHandleList(data)
            % constructor
            if nargin > 0
                obj.data = data;
            end
        end        
        
        function d = getData(obj,indices)
            % d = getData(obj,indices) returns the underlying data
            if nargin == 1
                indices = ':';
            end
            for o=1:numel(obj)
                d = obj(o).data(indices);
            end
        end
        
        function v = getInfo(obj,infofield)
            % v = getInfo(obj,infofield) returns information about the list
            % infofield is a optional parameter
            if nargin == 1
                v = [obj.info]; 
            else
                obj.checkForError('invalidInfoField',{infofield});
                v = arrayfun(@(o) o.info.(infofield),obj(:),'Uniform',0);
            end
        end
        
        function saveToFile(obj,filename)
            % saveToFile(obj,filename) saves the list to disk
            disp('Saving list, please wait...')
            dhl = obj.copy;
            save(filename,'dhl','-v7.3');
        end
        
        function n = getLength(obj)
            % n = getLength(obj) returns the size of the wrapped data
            n = cell2mat(arrayfun(@(o) numel(o.data),obj(:),'Uniform',0));
        end
        
        function s = getDataType(obj)
            % s = getDataType(obj) returns the class of the underlying 
            % data
            s = arrayfun(@(o) class(o.data),obj(:),'Uniform',0);
        end
        
        function dhl = flattenLists(obj)
            % dhl = flattenLists(obj) merges data of 2D list and returns a
            % new list. Metadata needs to be managed by the user
            merged_data = [];
            for o=1:numel(obj)
                merged_data = [merged_data; obj(o).getData]; %#ok<AGROW>
            end
            dhl = eval(strcat(class(obj),'(merged_data)')); 
        end
        
        function v = addFieldAndValues(obj,fieldname,values)
            % v = addFieldAndValues(obj,fieldname,values) adds a new field 
            % (and values for it) to the data struct.
            obj.checkForError('invalidObjectLength',{1});
            obj.checkForError('inputMismatch',{numel(values)});
            if ~iscell(values)
                values = num2cell(values);
            end
            [obj.data(:).(fieldname)] = deal(values{:});
            v = obj.get(fieldname);
        end
        
        function v = addFieldAndValuesById(obj,fieldname,values,idfield,ids)
            % v = addFieldAndValuesById(obj,fieldname,values,idfield,ids) 
            % adds a new field fieldname to the data struct and add values 
            % matched by id field will be assign values sorted by id
            % Example: 
            % dhl.addFieldAndValuesById('newfield',values,'names',ids)
            obj.checkForError('invalidObjectLength',{1});
            obj.checkForError('invalidFieldname',{idfield});
            obj.checkForError('duplicateFieldname',{fieldname});
            switch class(values)
                case 'cell'
                    d = repmat({''},obj.getLength,1);
                case 'double'
                    d = nan(obj.getLength,1);
            end
            [~,idx] = ismember(obj.get(idfield),ids);
            d(idx>0) = values(idx(idx>0));
            obj.addFieldAndValues(fieldname,d);
            v = obj.get(fieldname);
        end
        
        function removeField(obj,fieldname)
            % removeField(obj,fieldname) removes a field from the data 
            % struct
            for o=1:numel(obj)
                obj(o).checkForError('invalidFieldname',{fieldname});
                obj(o).data = rmfield(obj(o).data,fieldname);
            end
        end

        function v = get(obj,fieldname,indices)
            % get(obj,fieldname,indices) is used to get at the underlying
            % data of a GenericDataHandleList. Without parameters, it 
            % returns a list of matlab:hyperlinks that can be clicked to 
            % get to the respective fields. The first visible parameter is
            % the fieldname or the name of the method that will return a list
            % of results.The optional second parameter specifies indices if
            % only a subset of the data should be returned.
            if numel(obj) > 1
                warning('This method works only for 1x1 objects')
                return
            end
            if nargin == 1
                fns = fieldnames(obj.data);
                hist = com.mathworks.mlservices.MLCommandHistoryServices.getSessionHistory;
                lc = char(hist(end));
                for i=1:numel(fns)
                    cmd = strcat('<a href="matlab: ',lc(1:end-4),...
                        '.get(''',fns(i),''')">',fns(i),'</a>');
                    disp(cell2mat(cmd));
                end
                return
            elseif nargin == 2
                % numeric parameter -> return data for index
                if isnumeric(fieldname) || islogical(fieldname)
                    v = obj.getData(fieldname);
                    return
                else
                    indices = ':';
                end
            end
            if ismethod(obj,fieldname)
                v = obj.(fieldname)(indices);
            elseif isprop(obj,fieldname)
                v = obj.(fieldname);
            else
                fn =  strsplit(fieldname,'.');
                if length(fn) == 1
                    v = {obj.data(indices).(fieldname)}';
                else
                    v = obj.(cell2mat(fn(1))).(cell2mat(fn(2)));
                    v = {v(indices)'};
                end
            end
            if ~iscellstr(v) && iscell(v)
                if isobject(v{1})                    
                    v = [v{:}]';
                else
                    v = cell2mat(v);
                end
            end
        end

        function dhl = selectListsByInfoQuery(obj,varargin)
            % selectListsByInfoQuery(obj,varargin) creates a new array of
            % dhls selected from an array of source lists based on a query
            % Example: 
            % dhl.selectListsByInfoQuery({'name','=','Tom'})
            %
            % possible values for comparison:
            % equal:       '=','==','eq','is'
            % unequal:     '~','~=','!=','ne','not'
            % greater:     '>','gt','greater'; also '>=','ge'
            % lesser:      '<','lt','less'; also '<=','le'
            % between:     '><','between'
            % find in set: 'in','ismember'
            % find string: 'contains','like'
            % regexp:      'match','regexp'
            
            dhl = obj(obj.query(varargin{:},'info')).copy;
        end
        
        function l = query(obj,queryparams,level)
            % l = query(obj,queryparams,level) queries the data or the info
            % and returns a logial array
            if nargin == 2
                level = 'data';
            end
            n = round(numel(queryparams)/3);
            params = reshape(queryparams,3,n);
            keys = params(1,:);
            comp = params(2,:);
            vals = params(3,:);
            switch level
                case 'info'
                    idxmat = false(numel(obj),n);
                case 'data'
                    idxmat = false(obj.getLength,n);
            end
            for i=1:numel(comp)
                switch level
                case 'info'
                    res = obj.getInfo(keys{i});
                case 'data'
                    res = obj.get(keys{i});
                end
                if iscell(res)
                    if isnumeric(res{1})
                        res = cell2mat(res);
                    end
                end
                switch comp{i}
                    case {'=','==','eq','is'}
                        if isnumeric(res)
                            indices = (res == str2double(vals{i}));
                        elseif islogical(res{1})
                            if str2double(vals{i}) == 1
                                indices = cell2mat(res);
                            else
                                indices = ~cell2mat(res);
                            end
                        else
                            indices = strcmp(res,vals{i});
                        end
                    case {'~','~=','!=','ne','not'}
                        if isnumeric(res)
                            indices = (res ~= str2double(vals{i}));
                        elseif islogical(res{1})
                            if str2double(vals{i}) == 1
                                indices = ~cell2mat(res);
                            else
                                indices = cell2mat(res);
                            end
                        else
                            indices = ~strcmp(res,vals{i});    
                        end
                    case {'>','gt','greater'}
                        if isnumeric(res)
                            indices = (res > str2double(vals{i}));
                        else
                            warning([comp{i} ' ignored for non-numeric parameter ' keys{i}])
                        end
                    case {'>=','ge'}
                        if isnumeric(res)
                            indices = (res >= str2double(vals{i}));
                        else
                            warning([comp{i} ' ignored for non-numeric parameter ' keys{i}])
                        end
                    case {'<','lt','less'}
                        if isnumeric(res)
                            indices = (res < str2double(vals{i}));
                        else
                            warning([comp{i} ' ignored for non-numeric parameter ' keys{i}])
                        end
                    case {'<=','le'}
                        if isnumeric(res)
                            indices = (res <= str2double(vals{i}));
                        else
                            warning([comp{i} ' ignored for non-numeric parameter ' keys{i}])
                        end
                    case {'><','between'}
                        if isnumeric(res)
                            v = str2num(vals{i}); %#ok<ST2NM>
                            if numel(v) == 2                                
                                indices = (res > v(1)) & (res < v(2));
                            else
                                warning('Two values needed!')    
                            end
                        else
                            warning([comp{i} ' ignored for non-numeric parameter ' keys{i}])
                        end
                    case {'contains','like'}
                        indices = ~cellfun(@isempty,strfind(res,vals{i}));
                    case {'match','regexp'}
                        indices = ~cellfun(@isempty,regexp(res,vals{i}));
                    case {'in','ismember'}
                        if isnumeric(res)
                            indices = ismember(res,str2num(vals{i})); %#ok<ST2NM>
                        else
                            indices = ismember(res,vals{i});
                        end
                    otherwise
                        disp([comp{i} ' not recognized'])
                end
                idxmat(:,i) = indices;
            end
            l = logical(min(idxmat,[],2));
        end
        
        function indices = findDataByFieldValue(obj,fieldname,fieldvalue)
            % indices = findDataByFieldValue(obj,fieldname,fieldvalue) 
            % finds values, using hashed values.
            % If hash already exists for this fieldname, reuse
            obj.checkForError('invalidObjectLength',{1});
            if ~isfield(obj.hash,'fieldname')
                obj.buildHashForField(fieldname);
            elseif ~strcmp(obj.hash.fieldname,fieldname)
                obj.buildHashForField(fieldname);
            end
            [valid,key] = max(ismember(obj.hash.keys,fieldvalue));
            if valid == 0
                %key not found, set all indices to 0
                indices = false(numel(obj.hash.indices),1);
            else
                indices = ismember(obj.hash.indices,key);
            end
        end

        function removeDataByIndex(obj,indices)
            % removeDataByIndex(obj,indices) deletes data from underlying 
            % struct. Should be used with care
            obj.checkForError('invalidDataIndex',{indices});
            old_length = obj.getLength;
            obj.data(indices) = [];
            sprintf('old length: %d, new length: %d',old_length, obj.getLength)
        end
        
        function removeDataByFieldValue(obj,fieldname,fieldvalue)
            % removeDataByFieldValue(obj,fieldname,fieldvalue) deletes data 
            % from underlying struct, for data where specified field has a 
            % given value
            obj.checkForError('invalidObjectLength',{1});
            indices = obj.findDataByFieldValue(fieldname,fieldvalue);
            obj.removeDataByIndex(indices);
        end
        
        function sortDataByFieldValue(obj,fieldname,mode)
            % sortDataByFieldValue(obj,fieldname,mode) sorts the underlying
            % data by values for a field
            if nargin == 2
                mode = 'ascend';
            end            
            for o=1:numel(obj)
                [~, order] = sort(obj(o).get(fieldname),mode);
                obj(o).data = obj(o).data(order);
            end
        end
        
        function dhl = createSubsetByIndex(obj,indices) %#ok<INUSD>
            % dhl = createSubsetByIndex(obj,indices) creates a subset list 
            % from given indices
            dhl = eval(strcat(class(obj),'(obj.data(indices))'));
        end
        
        function dhls_array = divideListsByInfoField(obj,infofield)
            % dhls_array = divideListsByInfoField(obj,infofield) create a 
            % cell array of lists split by a values for an info field
            vals = obj.getInfo(infofield);
            keys = unique(vals);
            dhls_array = cell(numel(keys),1);
            for i = 1:length(keys)
                dhls_array{i} = obj(ismember(vals,keys(i)));
            end
        end
        
        function dhls = selectListsByInfoFieldValue(obj,infofield,value)
            % create a list selected by a info field value
            vals = obj.getInfo(infofield);
            idx = strcmp(vals,value);
            dhls = obj(idx);
        end

        function dhl = createSubsetByFieldValue(obj,fieldname,fieldvalue)
            % create a subset list (?dhl) for a given fieldvalue
            indices = obj.findDataByFieldValue(fieldname,fieldvalue);
            dhl = obj.createSubsetByIndex(indices);
        end
        
        function dhl = createSubsetByQuery(obj,queryparams,level)
            % createSubsetByQuery creates a subset list (?dhl) for a given 
            % query
            if nargin == 2
                level = 'data';
            end
            if numel(obj) > 1
                dhl = eval(strcat(class(obj),'.empty(',num2str(numel(obj)),',0)'));
            end
            for o=1:numel(obj)
                indices = obj(o).query(queryparams,level);
                dhl(o) = obj(o).createSubsetByIndex(indices);
            end
        end
       
        function dhls = splitByIndexMatrix(obj,matrix)
            % splitByIndexMatrix split list into multiple lists (?dhls) 
            % with a matrix of indices
            obj.checkForError('invalidObjectLength',{1});
            dhls = eval(strcat(class(obj),'.empty(size(indices,2),0);'));
            for i = 1:size(matrix,2)
                dhls(i) = obj.createSubsetByIndex(matrix(:,i));
            end
        end

        function dhls = splitByGroupIDs(obj,ids)
            % splitByGroupIDs split list into multiple lists (?dhls) by 
            % group IDs
            obj.checkForError('invalidObjectLength',{1});
            indices = obj.buildLookupTableForIDs(ids);
            dhls = eval(strcat(class(obj),'.empty(size(indices,2),0);'));
            for i = 1:size(indices,2)
                dhls(i) = obj.createSubsetByIndex(indices(:,i));
            end
        end
        
        function dhls = splitByFieldName(obj,fieldname)
            % splitByFieldName split list into multiple lists (?dhls) by 
            % fieldname
            obj.checkForError('invalidObjectLength',{1});
            indices = obj.buildLookupTableForField(fieldname);
            % use eval to handle base- and inherited classes
            dhls = eval(strcat(class(obj),'.empty(size(indices,2),0);'));
            for i = 1:size(indices,2)
                dhls(i) = obj.createSubsetByIndex(indices(:,i));
            end
        end
        
        function mergeInfoFields(obj,infofns,fnout)
            % creates strings!
            if nargin == 2                
                %fnout = strcat(infofns{:});
                fnout = strjoin(infofns,'_');
            end
            n = numel(infofns);
            cmat = cell(numel(obj),n);
            for i=1:n
                infodata = obj.getInfo(infofns{i});
                if isnumeric(infodata{1})
                    infodata = cellfun(@num2str,infodata,'un',0);
                end
                cmat(:,i) = infodata;
            end
            comVal = join(cmat,'-');%creates string!
            %convertStringsToChars 2018a
            for o=1:numel(obj)
                obj(o).info.(fnout) = comVal(o);
            end
        end        

        function addInfo(obj,fieldname,values)
            % addInfo(obj,fieldname,values) adds
            if numel(obj) == 1
                obj.info.(fieldname) = values;
            else
                 for o=1:numel(obj)
                     obj(o).info.(fieldname) = values{o};
                 end
            end
        end
        
        function addMetadataFromFile(obj,filename)
            % addMetadataFromFile(obj,filename) adds metadata saved in an
            % .xls file
            disp('Loading database, please wait...');
            obj.checkForError('fileNotFound',{filename});
            [~,sheets] = xlsfinfo(filename);
            db = [];
            for i=1:numel(sheets)
                sheet = cell2mat(sheets(i));
                [~,~,xlsdata] = xlsread(filename, sheet);
                db.(sheet) = xlsdata;
            end
            obj.info.dbpath = filename;
            %obj.applyDataBase;
            cns = fieldnames(db);
            for c=1:numel(cns)
                cname = cell2mat(cns(c));
                if ~isfield(obj.data,cname)
                    continue
                end
                ins = obj.db.(cname)(2:end,1);
                for i=1:numel(ins)
                    iname = cell2mat(ins(i));
                    fns = obj.db.(cname)(1,2:end);
                    indices = obj.findDataByFieldValue(cname,iname);
                    for f=1:numel(fns)
                        fname = cell2mat(fns(f));
                        b = cell2mat(obj.db.(cname)(i+1,f+1));
                        [obj.data(indices).(fname)] = deal(b);
                    end
                end
            end
        end
        
        function db = loadMetaData(obj,filename)
            disp('Loading database, please wait...');
            obj.checkForError('fileNotFound',{filename});
            [~,sheets] = xlsfinfo(filename);
            for i=1:numel(sheets)
                sheet = cell2mat(sheets(i));
                [~,~,xlsdata] = xlsread(filename, sheet);
                obj.db.(sheet) = xlsdata;
            end
            db = obj.db;
            obj.info.dbpath = filename;
        end
        
    end

    methods (Static)

        function dhl = LoadFromFile(filename)
            disp('Loading list from file, please wait...')
            if ~exist(filename,'file')
                error('GenericDataHandleList:fileNotFound',...
                    'Error. \nFile ''%s'' not found!',filename)
            end
            c = load(filename,'dhl');
            dhl = c.dhl;
            dhl.setup;
        end
        
    end
    
    methods (Access = protected)
        % hidden methods, for internal use only
        
        function this = findThis(obj) %#ok<MANU>
            hist = com.mathworks.mlservices.MLCommandHistoryServices.getSessionHistory;
            lastcom = char(hist(end));
            this = lastcom(1:strfind(lastcom,'.')-1);
        end
        
        function checkForError(obj,errorname,params)
            switch errorname
                case 'fileNotFound'
                    if ~exist(params{1},'file')
                        error('GenericDataHandleList:fileNotFound',...
                            'Error. \nFile ''%s'' not found!',params{1})
                    end
                case 'invalidObjectLength'
                    if numel(obj) ~= params{1}
                        error('GenericDataHandleList:invalidObjectLength',...
                            'Error. \nThis method can only be applied on 1x%d lists.',params{1})
                    end
                case 'differentObjectLength'
                    if numel(obj) ~= numel(params{1})
                        error('GenericDataHandleList:differentObjectLength',...
                            'Error. \nThe two objects differ in length.')
                    end
                case 'invalidDataIndex'
                    if islogical(params{1}) && (length(params{1}) ~= obj.getLength)
                        error('GenericDataHandleList:invalidDataIndex',...
                            'Error. \nWrong size of logical array.')
                    end
                    if isnumeric(params{1}) && (obj.getLength < max(params{1}))
                        error('GenericDataHandleList:invalidDataIndex',...
                            'Error. \nIndex exceeds length of data struct.')
                    end
                case 'nonscalarDataIndex'
                    if length(params{1}) > 1
                        error('GenericDataHandleList:nonscalarDataIndex',...
                            'Error. \nScalar index required. Please enter only a single number.')
                    end
                case 'invalidFieldname'
                    if ~isfield(obj.data,params{1})
                        error('GenericDataHandleList:invalidFieldname',...
                            'Error. \nA field named ''%s'' does not exist!',params{1})
                    end
                case 'duplicateFieldname'
                    if isfield(obj.data,params{1})
                        error('GenericDataHandleList:duplicateFieldname',...
                            'Error. \nA field named ''%s'' already exists!',params{1})
                    end
                case 'invalidInfoField'
                    if ~isfield(obj(1).info,params{1})%TODO: check size
                        error('GenericDataHandleList:invalidInfoField',...
                            'Error. \nA infofield named ''%s'' does not exist!',params{1})
                    end
                case 'invalidObjectProperty'
                    if ~isprop(obj,params{1})
                        error('GenericDataHandleList:invalidObjectProperty',...
                            'Error. \nA property named ''%s'' does not exist!',params{1})
                    end
                case 'differentListTypes'
                    if params{1} ~= 1
                        error('GenericDataHandleList:differentListTypes',...
                            'Error. \nDifferent List types found!')
                    end
                case 'inputMismatch'
                    if ~ismember(params{1}, [1 obj.getLength])
                        error('GenericDataHandleList:inputMismatch',...
                            'Error. \nNumber of values does not match!')
                    end
            end
        end
        
        function buildHashForField(obj,fieldname)
            v = obj.get(fieldname);
            if iscellstr(v) || ~iscell(v)
                [keys,values,indices] = unique(v);
            else
                [keys,values,indices] = unique(cell2mat(v));
            end
            obj.hash.keys = keys;
            obj.hash.values = values;
            obj.hash.indices = indices;
            obj.hash.fieldname = fieldname;
        end        
        
        function matrix = buildLookupTableForIDs(obj,ids) %#ok<INUSL>
            ymat = unique(ids);
            matrix = cast(zeros(length(ids),length(ymat)),'logical');
            for i = 1:length(ymat)
                matrix(:,i) = ismember(ids,ymat(i));
            end
        end
        
        function matrix = buildLookupTableForField(obj,fieldname)
            if ~isnumeric(obj.get(fieldname,1))
                obj.buildHashForField(fieldname);
                xmat = obj.hash.indices;
                ymat = 1:length(obj.hash.keys);
            else
                xmat = obj.get(fieldname);
                ymat = unique(obj.get(fieldname));
            end
            matrix = false(length(xmat),length(ymat));
            for i = 1:length(ymat)
                matrix(:,i) = ismember(xmat,ymat(i));
            end
        end
        
    end
    
end