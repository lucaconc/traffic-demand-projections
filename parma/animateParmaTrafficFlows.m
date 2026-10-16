function manifest = animateParmaTrafficFlows(R,varargin)
%ANIMATEPARMATRAFFICFLOWS Fast interactive Parma traffic-flow animator.
%
% Rebuilt from scratch for responsiveness. The GUI displays ONE Poisson slice
% at a time. All stored member paths contribute to the inferred directed-edge
% flow. Static road/flow graphics are rebuilt only when the selected slice
% changes. Each red marker is owned by one complete stored physical path and
% moves from that path's first node to its last node before wrapping to the
% beginning of the SAME path. Each timer tick updates only one marker object.
%
% Main options
%   DayType             '' | text label(s)
%   SeasonDate          'auto' | 'MM-DD' | datetime | 'all'
%   Year                'auto' | numeric | 'AllYears' | 'all'
%   PoissonSlice        'auto' | numeric stratum id/row
%   FrameRate           15 initial interactive FPS
%   MaxCars             250 nominal car budget
%   CarCountFactor      0.01 initial density fraction of MaxCars
%   VehicleSpeedKmH     120
%   VehicleSpeedFactor  10 visual speed multiplier, constant on all edges
%   FlowBins            4 number of batched edge-flow drawing groups
%   MergeDisplayEdges  true; merge compatible consecutive flow edges for drawing
%   MergeFlowTolerance  0.20 relative flow difference allowed while merging
%   MergeAngleToleranceDeg 15 maximum turn angle for one straight approximation
%   MergeMaxChainEdges  25 maximum original edges represented by one display segment
%   MaxDisplayedEdges   legacy compatibility only; all positive-flow edges are represented
%   ShowBackgroundRoads false; when true, draw merged display segments in gray underneath
%   ShowCars            true
%   ShowCarPaths        true; overlay the complete paths currently owning cars
%   Visible             'on'
%   FigureSize          [10 7]
%   PathWeightMode      'sampling' | 'equal'
%   Seed                1234
%
% Legacy compatibility options OutputFile/VideoFile are accepted but ignored;
% this fast interactive version intentionally does not export video.

ip=inputParser; ip.FunctionName=mfilename;
addParameter(ip,'DayType','',@isTextSelectionLocal);
addParameter(ip,'SeasonDate','auto',@(x)ischar(x)||isstring(x)||isdatetime(x));
addParameter(ip,'Year','auto',@(x)(ischar(x)||isstring(x))||(isnumeric(x)&&isscalar(x)&&isfinite(x)));
addParameter(ip,'PoissonSlice','auto',@(x)(ischar(x)||isstring(x))||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1));
addParameter(ip,'FrameRate',15,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x<=60);
addParameter(ip,'MaxInteractiveFrameRate',30,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x<=60);
addParameter(ip,'MaxCars',250,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==floor(x));
addParameter(ip,'CarCountFactor',0.01,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<=1);
addParameter(ip,'VehicleSpeedKmH',120,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'VehicleSpeedFactor',10,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'CarMarkerSize',6,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'CarPathSamples',96,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=16&&x<=512&&x==floor(x));
addParameter(ip,'FlowBins',4,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x<=20&&x==floor(x));
addParameter(ip,'MaxDisplayedEdges',80,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==floor(x)); % legacy, ignored
addParameter(ip,'MergeDisplayEdges',true,@isScalarLogicalLikeLocal);
addParameter(ip,'MergeFlowTolerance',0.20,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
addParameter(ip,'MergeAngleToleranceDeg',15,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=90);
addParameter(ip,'MergeMaxChainEdges',25,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==floor(x));
addParameter(ip,'ShowBackgroundRoads',false,@isScalarLogicalLikeLocal);
addParameter(ip,'FlowLineWidth',[0.6 5.0],@(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(1)>0&&x(2)>=x(1));
addParameter(ip,'RoadLineWidth',0.25,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'ShowCars',true,@isScalarLogicalLikeLocal);
addParameter(ip,'ShowCarPaths',true,@isScalarLogicalLikeLocal);
addParameter(ip,'Visible','on',@(x)ischar(x)||isstring(x));
addParameter(ip,'FigureSize',[10 7],@(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&all(x>0));
addParameter(ip,'PathWeightMode','sampling',@(x)ischar(x)||isstring(x));
addParameter(ip,'Seed',1234,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==floor(x));
addParameter(ip,'OutputFile','',@(x)ischar(x)||isstring(x));
addParameter(ip,'VideoFile','',@(x)ischar(x)||isstring(x));
% Accepted legacy options, intentionally unused in this fast fixed-slice GUI.
addParameter(ip,'DurationSeconds',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0));
addParameter(ip,'SecondsPerState',0.8,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'HoldSeconds',0.12,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(ip,'Loop',false,@isScalarLogicalLikeLocal);
addParameter(ip,'EnsureAllPositiveEdges',false,@isScalarLogicalLikeLocal);
addParameter(ip,'MinEdgeTraversalsPerSecond',0,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
parse(ip,varargin{:});
o=ip.Results;
% All output paths lie outside the library (see safeOutputPath).
o=guardOutputFields(o,{'OutputFile','VideoFile'});
o.DayType=string(o.DayType); o.DayType=o.DayType(:); o.DayType=strtrim(o.DayType); o.DayType=o.DayType(strlength(o.DayType)>0);
o.ShowCars=logical(o.ShowCars); o.ShowCarPaths=logical(o.ShowCarPaths);
o.ShowBackgroundRoads=logical(o.ShowBackgroundRoads);
o.MergeDisplayEdges=logical(o.MergeDisplayEdges);
o.PathWeightMode=lower(strtrim(char(string(o.PathWeightMode))));
if ~ismember(o.PathWeightMode,{'sampling','equal'})
    error('animateParmaTrafficFlows:PathWeightMode','PathWeightMode must be ''sampling'' or ''equal''.');
end
if strlength(string(o.OutputFile))>0 || strlength(string(o.VideoFile))>0
    warning('animateParmaTrafficFlows:NoVideo', ...
        'This fast interactive animator does not export video. OutputFile/VideoFile is ignored.');
end

[R,resultFile]=loadResultsLocal(R);
if ~isfield(R,'G') || (~isa(R.G,'graph') && ~isa(R.G,'digraph'))
    error('animateParmaTrafficFlows:Graph','Saved results must contain graph/digraph field R.G.');
end
if ~isfield(R,'muRatePerHourByStratum') || isempty(R.muRatePerHourByStratum)
    error('animateParmaTrafficFlows:Rates','R.muRatePerHourByStratum is required.');
end
if ~isfield(R,'strataTable') || isempty(R.strataTable)
    error('animateParmaTrafficFlows:Strata','R.strataTable is required.');
end
G=R.G; T=R.strataTable; M=orientRateMatrixLocal(R,T);
[T,M,selectedDayTypes]=filterDayTypesLocal(T,M,o.DayType);
[T,M,selectedYear]=filterYearLocal(T,M,o.Year);
[T,M,selectedSeasonDate,selectedSeasonDOY]=filterSeasonDateLocal(T,M,o.SeasonDate);
[T,M]=filterEstimatedStatesLocal(T,M);
[T,M]=sortAnimationStatesLocal(T,M);
M(~isfinite(M)|M<0)=0;
[p,nStates]=size(M);
if nStates<1, error('animateParmaTrafficFlows:Empty','No estimated Poisson slices remain after filtering.'); end
selectedSliceIndex=resolvePoissonSliceLocal(T,o.PoissonSlice);

[classPaths,pathSource,representativeFallback]=resolveClassPathsLocal(R,p);
for j=1:p
    classPaths{j}=cleanAndDeduplicatePathsLocal(classPaths{j},numedges(G));
end
[classPathWeights,pathWeightSource]=resolveClassPathWeightsLocal(R,classPaths,p,o.PathWeightMode);
classIDs=find(cellfun(@(x)~isempty(x),classPaths));
if isempty(classIDs), error('animateParmaTrafficFlows:NoPaths','No drawable member paths are stored in the results.'); end

% Flatten the retained physical member paths once.  Cars are allocated to
% THESE full paths (not to edges or merged display segments).  Path q carries
% M(class(q),slice)*weight(q) vehicles/hour in the selected Poisson slice.
[flatPaths,flatPathClass,flatPathWeight]=flattenPathCatalogLocal(classPaths,classPathWeights,classIDs);
if isempty(flatPaths), error('animateParmaTrafficFlows:NoPaths','No physical member paths remain after cleaning.'); end

% Exact directed-edge flow for every slice. This is the only expensive path
% aggregation and is done once before the GUI starts.
edgeFlowAll=buildEdgeFlowMatrixLocal(classPaths,classPathWeights,classIDs,M,numedges(G));
activeEdgeIDs=find(max(edgeFlowAll,[],2)>0);
edgeFlowAll=edgeFlowAll(activeEdgeIDs,:);
edgeGeom=directedEdgeGeometryLocal(G,activeEdgeIDs);
valid=isfinite(edgeGeom.lengthKm)&edgeGeom.lengthKm>0&isfinite(edgeGeom.lat1)&isfinite(edgeGeom.lon1)&isfinite(edgeGeom.lat2)&isfinite(edgeGeom.lon2);
activeEdgeIDs=activeEdgeIDs(valid); edgeFlowAll=edgeFlowAll(valid,:); edgeGeom=subsetEdgeGeometryLocal(edgeGeom,valid);
if isempty(activeEdgeIDs), error('animateParmaTrafficFlows:NoEdges','No positive-flow directed edges are drawable.'); end

% GUI ---------------------------------------------------------------------
f=figure('Name','Parma traffic flow - fast animator','Color','w','Visible',char(string(o.Visible)), ...
    'Units','inches','Position',[1 1 double(o.FigureSize(:).')]);
try, set(f,'GraphicsSmoothing','off'); catch, end
ax=axes(f,'Units','normalized','Position',[0.06 0.20 0.91 0.75]); hold(ax,'on');
% The full positive-flow network is represented, but compatible consecutive
% edges are merged into longer display segments when each slice is selected.
% This keeps topology/flow coverage while reducing rendered vertex count.
roadHandle=plot(ax,nan,nan,'Color',[0.86 0.86 0.86], ...
    'LineWidth',o.RoadLineWidth,'HandleVisibility','off');
axis(ax,'equal'); box(ax,'on'); grid(ax,'off');
lonAll=[edgeGeom.lon1(:);edgeGeom.lon2(:)]; latAll=[edgeGeom.lat1(:);edgeGeom.lat2(:)];
lonMin=min(lonAll); lonMax=max(lonAll); latMin=min(latAll); latMax=max(latAll);
dLon=max(lonMax-lonMin,1e-4); dLat=max(latMax-latMin,1e-4);
xlim(ax,[lonMin-0.03*dLon lonMax+0.03*dLon]); ylim(ax,[latMin-0.03*dLat latMax+0.03*dLat]);
xlabel(ax,'Longitude'); ylabel(ax,'Latitude');
set(ax,'FontName','Helvetica','FontSize',9,'TickDir','out');

% Fixed number of batched flow objects. Only their X/Y data changes on slice
% selection; they are never touched by the animation timer.
flowHandles=gobjects(o.FlowBins,1);
flowColors=parula(max(o.FlowBins,2));
for b=1:o.FlowBins
    flowHandles(b)=plot(ax,nan,nan,'Color',flowColors(b,:), ...
        'LineWidth',o.FlowLineWidth(1)+(o.FlowLineWidth(2)-o.FlowLineWidth(1))*(b-1)/max(o.FlowBins-1,1), ...
        'HandleVisibility','off');
end

% One batched overlay shows the COMPLETE physical paths that own visible cars.
% It is rebuilt only when slice/density changes and is never touched by timer.
carPathHandle=plot(ax,nan,nan,'LineStyle','-','Color',[0.85 0.45 0.45], ...
    'LineWidth',0.75,'HandleVisibility','off');
if ~o.ShowCarPaths, set(carPathHandle,'Visible','off'); end

% One lightweight car graphics object. No halo, no per-car graphics objects.
carHandle=plot(ax,nan,nan,'LineStyle','none','Marker','o','MarkerSize',o.CarMarkerSize, ...
    'MarkerFaceColor',[0.90 0.10 0.02],'MarkerEdgeColor','none','HandleVisibility','off');
if ~o.ShowCars, set(carHandle,'Visible','off'); end
statusText=text(ax,0.012,0.018,'','Units','normalized','VerticalAlignment','bottom', ...
    'FontSize',8,'BackgroundColor','w','Margin',2,'Clipping','on');
mainTitle=title(ax,'');

sliceLabels=poissonSliceLabelsLocal(T);
uicontrol(f,'Style','text','Units','normalized','Position',[0.06 0.115 0.105 0.035], ...
    'String','Poisson slice','BackgroundColor','w','HorizontalAlignment','left','FontWeight','bold');
slicePopup=uicontrol(f,'Style','popupmenu','Units','normalized','Position',[0.165 0.112 0.50 0.043], ...
    'String',cellstr(sliceLabels),'Value',selectedSliceIndex,'Callback',@onSliceChanged);
playButton=uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.69 0.112 0.10 0.043], ...
    'String','Pause','FontWeight','bold','Callback',@onPlayPause);
stopButton=uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.80 0.112 0.10 0.043], ...
    'String','Stop','FontWeight','bold','Callback',@onStop);

carsEnabled=o.ShowCars;
carCountFactorLive=double(o.CarCountFactor);
frameRateLive=min(double(o.FrameRate),double(o.MaxInteractiveFrameRate));
uicontrol(f,'Style','text','Units','normalized','Position',[0.06 0.055 0.06 0.03], ...
    'String','Cars','BackgroundColor','w','HorizontalAlignment','left','FontWeight','bold');
carsPopup=uicontrol(f,'Style','popupmenu','Units','normalized','Position',[0.12 0.052 0.12 0.038], ...
    'String',{'Cars On','Cars Off'},'Value',1+double(~carsEnabled),'Callback',@onCarsChanged);
uicontrol(f,'Style','text','Units','normalized','Position',[0.27 0.055 0.075 0.03], ...
    'String','Density','BackgroundColor','w','HorizontalAlignment','left','FontWeight','bold');
densitySlider=uicontrol(f,'Style','slider','Units','normalized','Position',[0.345 0.058 0.22 0.026], ...
    'Min',0,'Max',1,'Value',densityFactorToSliderLocal(carCountFactorLive,1e-4,1), ...
    'Callback',@onDensityChanged);
densityText=uicontrol(f,'Style','text','Units','normalized','Position',[0.57 0.055 0.09 0.03], ...
    'String',sprintf('%.2f%%',100*carCountFactorLive),'BackgroundColor','w','HorizontalAlignment','left');
uicontrol(f,'Style','text','Units','normalized','Position',[0.68 0.055 0.05 0.03], ...
    'String','FPS','BackgroundColor','w','HorizontalAlignment','left','FontWeight','bold');
fpsSlider=uicontrol(f,'Style','slider','Units','normalized','Position',[0.73 0.058 0.17 0.026], ...
    'Min',5,'Max',30,'Value',min(30,max(5,frameRateLive)),'Callback',@onFPSChanged);
fpsText=uicontrol(f,'Style','text','Units','normalized','Position',[0.905 0.055 0.05 0.03], ...
    'String',sprintf('%d',round(frameRateLive)),'BackgroundColor','w','HorizontalAlignment','left');

% Dynamic state: every visible car is permanently associated with one FULL
% physical member path until the slice/density selection is rebuilt.  The
% timer changes only carPhase and the marker X/Y coordinates.
carPathRow=zeros(0,1); carPhysicalPathID=zeros(0,1); carPhase=zeros(0,1); carInvLength=zeros(0,1);
carLatLookup=zeros(0,double(o.CarPathSamples)); carLonLookup=zeros(0,double(o.CarPathSamples));
displayGeom=emptyDisplayGeomLocal(); displayFlow=zeros(0,1); displayCarWeight=zeros(0,1);
effectiveSpeed=double(o.VehicleSpeedKmH)*double(o.VehicleSpeedFactor);
lastTick=tic; stopped=false;
rebuildSlice();

timerObj=timer('ExecutionMode','fixedSpacing','Period',1/frameRateLive,'BusyMode','drop', ...
    'TimerFcn',@onTick,'ErrorFcn',@onTimerError);
set(f,'CloseRequestFcn',@onClose);
start(timerObj);

manifest=struct();
manifest.resultFile=resultFile;
manifest.figure=f; manifest.timer=timerObj;
manifest.edgeFlowByState=edgeFlowAll; manifest.activeEdgeIDs=activeEdgeIDs;
manifest.pathSource=pathSource; manifest.pathWeightSource=pathWeightSource;
manifest.carMotion='full-physical-path'; manifest.numberOfPhysicalPaths=numel(flatPaths); manifest.showCarPaths=o.ShowCarPaths;
manifest.maxDisplayedEdges=o.MaxDisplayedEdges; manifest.showBackgroundRoads=o.ShowBackgroundRoads;
manifest.mergeDisplayEdges=o.MergeDisplayEdges; manifest.mergeFlowTolerance=o.MergeFlowTolerance;
manifest.mergeAngleToleranceDeg=o.MergeAngleToleranceDeg; manifest.mergeMaxChainEdges=o.MergeMaxChainEdges; manifest.carPathSamples=o.CarPathSamples;
manifest.dayTypes=selectedDayTypes; manifest.year=selectedYear;
manifest.seasonDate=selectedSeasonDate; manifest.seasonDateDOY=selectedSeasonDOY;
manifest.controls=struct('slicePopup',slicePopup,'carsPopup',carsPopup,'densitySlider',densitySlider,'fpsSlider',fpsSlider);
if any(representativeFallback(classIDs))
    warning('animateParmaTrafficFlows:RepresentativeFallback', ...
        'Some classes only contain representative paths; unavailable alternatives cannot be visualized.');
end

    function rebuildSlice()
        selectedSliceIndex=get(slicePopup,'Value');
        rate=edgeFlowAll(:,selectedSliceIndex);
        rate(~isfinite(rate)|rate<0)=0;

        % Every positive-flow original edge is represented.  For rendering
        % speed, compatible 1-in/1-out chains are approximated by one straight
        % display segment.  The original edge-flow matrix is never modified.
        if o.MergeDisplayEdges
            [displayGeom,displayFlow,displayCarWeight,mergeInfo]=mergeDisplayEdgesLocal( ...
                edgeGeom,rate,double(o.MergeFlowTolerance),double(o.MergeAngleToleranceDeg),double(o.MergeMaxChainEdges));
        else
            pos=find(rate>0);
            displayGeom=displayGeomFromOriginalLocal(edgeGeom,pos);
            displayFlow=rate(pos);
            displayCarWeight=rate(pos);
            mergeInfo=struct('originalPositiveEdges',numel(pos),'displaySegments',numel(pos),'mergedEdgesSaved',0);
        end
        % Merged edge geometry is ONLY a rendering optimization for the static
        % flow layer.  Vehicle trajectories are built separately from the full
        % stored physical paths below.
        updateFlowBins(displayGeom,displayFlow);
        if o.ShowBackgroundRoads && ~isempty(displayFlow)
            [rLat,rLon]=displaySegmentsLatLonLocal(displayGeom,1:numel(displayFlow));
            set(roadHandle,'XData',rLon,'YData',rLat);
        else
            set(roadHandle,'XData',nan,'YData',nan);
        end

        % Allocate cars to COMPLETE physical member paths.  The flow of a
        % physical path is its class Poisson rate times its within-class path
        % weight.  Therefore every red dot has an explicit path owner.
        pathRate=double(M(flatPathClass,selectedSliceIndex)).*flatPathWeight;
        pathRate(~isfinite(pathRate)|pathRate<0)=0;
        nTarget=max(1,round(double(o.MaxCars)*carCountFactorLive));
        counts=proportionalIntegerCountsLocal(pathRate,nTarget,false);
        ownedPathIDs=find(counts>0);
        [ownedLatLUT,ownedLonLUT,ownedLengthKm,ownedValid]=samplePhysicalPathsLocal( ...
            G,flatPaths(ownedPathIDs),double(o.CarPathSamples));
        % Drop any pathological stored path whose ordered edge geometry could
        % not be decoded.  This should normally be empty.
        if any(~ownedValid)
            badIDs=ownedPathIDs(~ownedValid); counts(badIDs)=0; ownedPathIDs=ownedPathIDs(ownedValid);
            ownedLatLUT=ownedLatLUT(ownedValid,:); ownedLonLUT=ownedLonLUT(ownedValid,:); ownedLengthKm=ownedLengthKm(ownedValid);
        end
        n=sum(counts);
        carPathRow=zeros(n,1); carPhysicalPathID=zeros(n,1); carPhase=zeros(n,1); carInvLength=zeros(n,1);
        k=0;
        for rr=1:numel(ownedPathIDs)
            pid=ownedPathIDs(rr); c=counts(pid); if c<=0, continue; end
            ids=k+(1:c); k=k+c;
            carPathRow(ids)=rr; carPhysicalPathID(ids)=pid;
            % The first dot on every path starts exactly at the path origin.
            % Additional dots on the same path are evenly spaced so they remain
            % separately visible, but every one wraps through origin -> end.
            carPhase(ids)=(0:c-1).'/c;
            carInvLength(ids)=1/max(ownedLengthKm(rr),eps);
        end
        carLatLookup=ownedLatLUT; carLonLookup=ownedLonLUT;
        if o.ShowCarPaths && ~isempty(carLonLookup)
            [cpLat,cpLon]=lookupRowsToLineLocal(carLatLookup,carLonLookup);
            set(carPathHandle,'XData',cpLon,'YData',cpLat,'Visible','on');
        else
            set(carPathHandle,'XData',nan,'YData',nan);
        end
        updateCarXY();
        if carsEnabled, set(carHandle,'Visible','on'); else, set(carHandle,'Visible','off'); end
        label=stratumLabelLocal(T,selectedSliceIndex);
        set(mainTitle,'String',sprintf('Poisson slice %s | %s',poissonSliceIdLocal(T,selectedSliceIndex),label));
        set(statusText,'String',sprintf('%d positive-flow edges -> %d merged display segments | %d cars owned by %d full paths | speed: %.0f km/h | FPS: %d', ...
            mergeInfo.originalPositiveEdges,mergeInfo.displaySegments,n,numel(ownedPathIDs),effectiveSpeed,round(frameRateLive)));
        drawnow limitrate;
    end

    function updateFlowBins(g,rates)
        for bb=1:o.FlowBins, set(flowHandles(bb),'XData',nan,'YData',nan); end
        if isempty(rates), return; end
        rates=double(rates(:)); pos=rates(rates>0);
        if isempty(pos), return; end
        lo=min(pos); hi=max(pos);
        if hi<=lo
            localBin=ones(size(rates)); localBin(rates<=0)=0;
        else
            z=(log1p(rates)-log1p(lo))/max(log1p(hi)-log1p(lo),eps);
            localBin=1+floor(z*o.FlowBins); localBin(rates<=0)=0; localBin=min(localBin,o.FlowBins);
        end
        for bb=1:o.FlowBins
            ids=find(localBin==bb);
            if isempty(ids), continue; end
            [lat,lon]=displaySegmentsLatLonLocal(g,ids);
            set(flowHandles(bb),'XData',lon,'YData',lat);
        end
    end

    function updateCarXY()
        if isempty(carPathRow) || isempty(carLonLookup)
            set(carHandle,'XData',nan,'YData',nan); return;
        end
        % Each row in the lookup tables is one COMPLETE stored physical path.
        % carPathRow never changes during playback, so a red dot cannot jump
        % from one route to another.
        e=carPathRow; u=min(max(carPhase,0),1-eps);
        nS=size(carLonLookup,2);
        z=u*(nS-1); j0=floor(z)+1; a=z-floor(z); j1=min(j0+1,nS);
        ind0=sub2ind(size(carLonLookup),e,j0);
        ind1=sub2ind(size(carLonLookup),e,j1);
        lon=(1-a).*carLonLookup(ind0)+a.*carLonLookup(ind1);
        lat=(1-a).*carLatLookup(ind0)+a.*carLatLookup(ind1);
        good=isfinite(lon)&isfinite(lat); lon(~good)=nan; lat(~good)=nan;
        set(carHandle,'XData',lon,'YData',lat);
    end

    function onTick(~,~)
        if ~isgraphics(f) || stopped || ~carsEnabled || isempty(carPathRow), lastTick=tic; return; end
        dt=toc(lastTick); lastTick=tic;
        if ~(isfinite(dt)&&dt>0), dt=1/frameRateLive; end
        % Same physical distance for every car; long edges therefore take
        % proportionally longer than short edges.
        carPhase=mod(carPhase + (effectiveSpeed*dt/3600).*carInvLength,1);
        updateCarXY();
        drawnow limitrate nocallbacks;
    end

    function onSliceChanged(~,~)
        stopped=false; set(playButton,'String','Pause');
        rebuildSlice(); lastTick=tic;
        if strcmp(timerObj.Running,'off'), start(timerObj); end
    end

    function onCarsChanged(src,~)
        carsEnabled=(get(src,'Value')==1);
        if carsEnabled
            set(carHandle,'Visible','on'); stopped=false; set(playButton,'String','Pause'); lastTick=tic;
            if strcmp(timerObj.Running,'off'), start(timerObj); end
        else
            set(carHandle,'Visible','off');
        end
    end

    function onDensityChanged(src,~)
        carCountFactorLive=sliderToDensityFactorLocal(get(src,'Value'),1e-4,1);
        set(densityText,'String',sprintf('%.2f%%',100*carCountFactorLive));
        rebuildSlice(); lastTick=tic;
    end

    function onFPSChanged(src,~)
        newFPS=round(get(src,'Value')); newFPS=max(5,min(30,newFPS));
        frameRateLive=newFPS; set(fpsText,'String',sprintf('%d',newFPS));
        wasRunning=strcmp(timerObj.Running,'on'); if wasRunning, stop(timerObj); end
        timerObj.Period=1/newFPS; lastTick=tic;
        if wasRunning && ~stopped, start(timerObj); end
        set(statusText,'String',sprintf('FPS changed to %d. Cars update only; static map is not redrawn.',newFPS));
    end

    function onPlayPause(src,~)
        if strcmp(timerObj.Running,'on')
            stop(timerObj); set(src,'String','Play');
        else
            stopped=false; lastTick=tic; start(timerObj); set(src,'String','Pause');
        end
    end

    function onStop(~,~)
        if strcmp(timerObj.Running,'on'), stop(timerObj); end
        stopped=true;
        % Put the first marker of each owned path at its path origin and keep
        % any additional markers evenly spaced on that same path.
        carPhase=resetOwnedPathPhasesLocal(carPathRow);
        updateCarXY(); drawnow limitrate; set(playButton,'String','Play');
    end

    function onTimerError(~,evt)
        try, warning('animateParmaTrafficFlows:Timer','Timer error: %s',evt.Data.Message); catch, end
    end

    function onClose(~,~)
        try, if isvalid(timerObj), stop(timerObj); delete(timerObj); end, catch, end
        delete(f);
    end
end

% -------------------------------------------------------------------------
% Physical-path car helpers.
function [flatPaths,flatClass,flatWeight]=flattenPathCatalogLocal(classPaths,classPathWeights,classIDs)
flatPaths=cell(0,1); flatClass=zeros(0,1); flatWeight=zeros(0,1);
for ii=1:numel(classIDs)
    j=classIDs(ii); pp=classPaths{j}; ww=classPathWeights{j};
    for k=1:numel(pp)
        if isempty(pp{k}), continue; end
        flatPaths{end+1,1}=double(pp{k}(:).'); %#ok<AGROW>
        flatClass(end+1,1)=j; %#ok<AGROW>
        if numel(ww)>=k && isfinite(ww(k)) && ww(k)>0, flatWeight(end+1,1)=double(ww(k)); else, flatWeight(end+1,1)=1/max(numel(pp),1); end %#ok<AGROW>
    end
end
end

function [latLUT,lonLUT,lengthKm,valid]=samplePhysicalPathsLocal(G,paths,nSamples)
% Convert each ordered edge-ID path into a cumulative-distance polyline and
% resample it at equal distance fractions.  A row is one immutable car route.
n=numel(paths); nSamples=max(2,round(nSamples));
latLUT=nan(n,nSamples); lonLUT=nan(n,nSamples); lengthKm=nan(n,1); valid=false(n,1);
phase=linspace(0,1,nSamples);
nodeLat=double(G.Nodes.Lat); nodeLon=double(G.Nodes.Lon);
for q=1:n
    e=round(double(paths{q}(:))); e=e(isfinite(e)&e>=1&e<=numedges(G));
    if isempty(e), continue; end
    [s,t]=edgeEndpointIndicesLocal(G,e);
    % Preserve the saved edge order. Normally t(k)==s(k+1). If a stored path
    % has a small discontinuity, bridge it explicitly rather than dropping the
    % rest of the path, so the marker still progresses from first to last edge.
    latv=nodeLat(s(1)); lonv=nodeLon(s(1));
    prev=t(1); latv(end+1,1)=nodeLat(prev); lonv(end+1,1)=nodeLon(prev); %#ok<AGROW>
    for k=2:numel(e)
        if s(k)~=prev
            latv(end+1,1)=nodeLat(s(k)); lonv(end+1,1)=nodeLon(s(k)); %#ok<AGROW>
        end
        latv(end+1,1)=nodeLat(t(k)); lonv(end+1,1)=nodeLon(t(k)); %#ok<AGROW>
        prev=t(k);
    end
    seg=haversineVectorKmLocal(latv(1:end-1),lonv(1:end-1),latv(2:end),lonv(2:end));
    keepSeg=isfinite(seg)&seg>1e-9;
    if ~any(keepSeg), continue; end
    % Remove zero-length vertices while preserving order.
    keep=[true; keepSeg(:)]; latv=latv(keep); lonv=lonv(keep);
    if numel(latv)<2, continue; end
    seg=haversineVectorKmLocal(latv(1:end-1),lonv(1:end-1),latv(2:end),lonv(2:end));
    cum=[0;cumsum(seg(:))]; L=cum(end);
    if ~(isfinite(L)&&L>0), continue; end
    target=phase*L;
    latLUT(q,:)=interp1(cum,latv,target,'linear');
    lonLUT(q,:)=interp1(cum,lonv,target,'linear');
    lengthKm(q)=L; valid(q)=all(isfinite(latLUT(q,:)))&&all(isfinite(lonLUT(q,:)));
end
end

function phase=resetOwnedPathPhasesLocal(carPathRow)
phase=zeros(size(carPathRow));
if isempty(carPathRow), return; end
u=unique(carPathRow(:).','stable');
for r=u
    ids=find(carPathRow==r); c=numel(ids); phase(ids)=(0:c-1).'/max(c,1);
end
end

function [latLine,lonLine]=lookupRowsToLineLocal(latLUT,lonLUT)
% Convert path lookup rows into one NaN-separated polyline object.
if isempty(latLUT), latLine=zeros(0,1); lonLine=zeros(0,1); return; end
[n,m]=size(latLUT); latLine=nan(n*(m+1),1); lonLine=nan(n*(m+1),1);
for q=1:n
    base=(q-1)*(m+1); latLine(base+(1:m))=latLUT(q,:); lonLine(base+(1:m))=lonLUT(q,:);
end
end

% -------------------------------------------------------------------------
% Helpers retained only for result compatibility and path/stratum decoding.
function v=densityFactorToSliderLocal(f,minF,maxF)
f=min(max(double(f),minF),maxF); v=(log10(f)-log10(minF))/(log10(maxF)-log10(minF));
end
function f=sliderToDensityFactorLocal(v,minF,maxF)
v=min(max(double(v),0),1); f=10^(log10(minF)+v*(log10(maxF)-log10(minF)));
end
function [T,M,selected]=filterDayTypesLocal(T,M,wanted)
if isempty(wanted) || any(strcmpi(wanted,'all')), selected="all"; return; end
if ~ismember('DayType',T.Properties.VariableNames), error('animateParmaTrafficFlows:DayType','strataTable has no DayType field.'); end
x=string(T.DayType); keep=false(height(T),1); for k=1:numel(wanted), keep=keep|strcmpi(x,wanted(k)); end
if ~any(keep), error('animateParmaTrafficFlows:DayType','Requested DayType not found. Available: %s',strjoin(unique(x),', ')); end
T=T(keep,:); M=M(:,keep); selected=strjoin(unique(x(keep),'stable'),', ');
end
function [T,M,selected]=filterYearLocal(T,M,wanted)
if ischar(wanted)||isstring(wanted)
    s=string(wanted); if strcmpi(s,'all'), selected="all"; return; end
    if strcmpi(s,'auto')
        if ismember('Year',T.Properties.VariableNames)
            y=double(T.Year); fy=y(isfinite(y)&y>0); if ~isempty(fy), wanted=max(fy); else, selected="auto"; return; end
        else, selected="auto"; return; end
    elseif strcmpi(s,'AllYears')
        if ismember('YearLabel',T.Properties.VariableNames), keep=strcmpi(string(T.YearLabel),'AllYears'); else, keep=true(height(T),1); end
        if ~any(keep), error('animateParmaTrafficFlows:Year','AllYears not found.'); end
        T=T(keep,:); M=M(:,keep); selected="AllYears"; return
    else
        z=str2double(s); if isfinite(z), wanted=z; else, error('animateParmaTrafficFlows:Year','Unknown Year selection %s.',s); end
    end
end
if isnumeric(wanted)
    if ~ismember('Year',T.Properties.VariableNames), error('animateParmaTrafficFlows:Year','strataTable has no Year field.'); end
    keep=double(T.Year)==double(wanted); if ~any(keep), error('animateParmaTrafficFlows:Year','Year %g not found.',wanted); end
    T=T(keep,:); M=M(:,keep); selected=string(wanted);
end
end
function [T,M,label,doy]=filterSeasonDateLocal(T,M,wanted)
label="all"; doy=NaN;
if ischar(wanted)||isstring(wanted)
    s=string(wanted); if strcmpi(s,'all'), return; end
    if strcmpi(s,'auto')
        if ismember('SeasonStartDOY',T.Properties.VariableNames)&&ismember('SeasonEndDOY',T.Properties.VariableNames)
            candidates=unique([double(T.SeasonStartDOY); double(T.SeasonEndDOY)]); candidates=candidates(isfinite(candidates));
            if isempty(candidates), return; end
            best=candidates(1); bestN=-1;
            for d=candidates(:).', k=seasonContainsLocal(double(T.SeasonStartDOY),double(T.SeasonEndDOY),d); if nnz(k)>bestN, best=d; bestN=nnz(k); end, end
            doy=best; label=monthDayLabelLocal(doy);
        else, return; end
    else
        doy=monthDayToDOYLocal(s); label=s;
    end
elseif isdatetime(wanted)
    doy=day(wanted(1),'dayofyear'); label=string(datestr(wanted(1),'mm-dd'));
end
if isfinite(doy)&&ismember('SeasonStartDOY',T.Properties.VariableNames)&&ismember('SeasonEndDOY',T.Properties.VariableNames)
    keep=seasonContainsLocal(double(T.SeasonStartDOY),double(T.SeasonEndDOY),doy);
    if ~any(keep), error('animateParmaTrafficFlows:SeasonDate','No Poisson slice contains SeasonDate %s.',label); end
    T=T(keep,:); M=M(:,keep);
end
end
function keep=seasonContainsLocal(a,b,d)
keep=(a<=b & d>=a & d<=b) | (a>b & (d>=a | d<=b));
end
function doy=monthDayToDOYLocal(v)
s=char(string(v)); parts=sscanf(s,'%d-%d'); if numel(parts)~=2, error('animateParmaTrafficFlows:SeasonDate','SeasonDate must be MM-DD, auto, all, or datetime.'); end
d=datetime(2001,parts(1),parts(2)); doy=day(d,'dayofyear');
end
function s=monthDayLabelLocal(doy)
d=datetime(2001,1,1)+days(round(doy)-1); s=string(datestr(d,'mm-dd'));
end
function [T,M]=filterEstimatedStatesLocal(T,M)
if ismember('Estimated',T.Properties.VariableNames), keep=logical(T.Estimated); else, keep=all(isfinite(M),1).'; end
T=T(keep,:); M=M(:,keep);
end
function [T,M]=sortAnimationStatesLocal(T,M)
vars=T.Properties.VariableNames; keys={};
for v={'StartMinute','EndMinute','Stratum'}, if ismember(v{1},vars), keys{end+1}=v{1}; end, end %#ok<AGROW>
if ~isempty(keys), [T,ord]=sortrows(T,keys); M=M(:,ord); end
end
function M=orientRateMatrixLocal(R,T)
M=double(R.muRatePerHourByStratum); n=height(T);
if size(M,2)==n, return; elseif size(M,1)==n, M=M.'; else, error('animateParmaTrafficFlows:Rates','Rate matrix dimensions do not match strataTable.'); end
end
function [classPaths,source,fallbackMask]=resolveClassPathsLocal(R,p)
classPaths=cell(p,1); fallbackMask=false(p,1); source='none';
fields={'classEdgePaths','memberEdgePaths','measurementClassEdgePaths','physicalPathsByClass'};
for q=1:numel(fields)
    if isfield(R,fields{q})&&~isempty(R.(fields{q})), classPaths=normalizeClassIndexedPathsLocal(R.(fields{q}),p); source=fields{q}; break; end
end
if all(cellfun(@isempty,classPaths)) && isfield(R,'representativeEdgePaths')
    v=R.representativeEdgePaths; if iscell(v), for j=1:min(p,numel(v)), if ~isempty(v{j}), classPaths{j}={double(v{j}(:).')}; fallbackMask(j)=true; end, end, source='representativeEdgePaths'; end
end
% Fill any missing class with representative fallback when available.
if isfield(R,'representativeEdgePaths')&&iscell(R.representativeEdgePaths)
    for j=1:min(p,numel(R.representativeEdgePaths))
        if isempty(classPaths{j})&&~isempty(R.representativeEdgePaths{j}), classPaths{j}={double(R.representativeEdgePaths{j}(:).')}; fallbackMask(j)=true; end
    end
end
end
function out=normalizeClassIndexedPathsLocal(v,p)
out=cell(p,1); if ~iscell(v), return; end
for j=1:min(p,numel(v)), out{j}=normalizeOneClassPathsLocal(v{j}); end
end
function paths=normalizeOneClassPathsLocal(v)
paths={}; if isempty(v), return; end
if isnumeric(v), if isvector(v), paths={double(v(:).')}; else, for r=1:size(v,1), z=v(r,:); z=z(isfinite(z)&z>0); if ~isempty(z), paths{end+1}=double(z(:).'); end, end, end
elseif iscell(v), for k=1:numel(v), if isnumeric(v{k})&&~isempty(v{k}), paths{end+1}=double(v{k}(:).'); end, end
end
end
function paths=cleanAndDeduplicatePathsLocal(paths,nEdges)
out={}; keys=strings(0,1);
for k=1:numel(paths), e=round(double(paths{k}(:).')); e=e(isfinite(e)&e>=1&e<=nEdges); if isempty(e), continue; end
    key=strjoin(string(e),','); if ~any(keys==key), out{end+1}=e; keys(end+1)=key; end %#ok<AGROW>
end
paths=out;
end
function [weights,source]=resolveClassPathWeightsLocal(R,classPaths,p,mode)
weights=cell(p,1); useSampling=strcmpi(mode,'sampling'); hitFields={'classPathHits','memberPathHits','classMemberPathHits'}; hits=[];
if useSampling
    for k=1:numel(hitFields), if isfield(R,hitFields{k})&&~isempty(R.(hitFields{k})), hits=R.(hitFields{k}); source=hitFields{k}; break; end, end
else, source='equal'; end
if ~exist('source','var'), source='equal-fallback'; end
for j=1:p
    n=numel(classPaths{j}); if n==0, weights{j}=zeros(0,1); continue; end
    w=ones(n,1);
    if useSampling&&~isempty(hits)&&iscell(hits)&&j<=numel(hits)&&~isempty(hits{j})
        h=double(hits{j}(:)); if numel(h)>=n, h=h(1:n); h(~isfinite(h)|h<0)=0; if sum(h)>0, w=h; end, end
    end
    weights{j}=w/sum(w);
end
end
function edgeFlow=buildEdgeFlowMatrixLocal(classPaths,classPathWeights,classIDs,M,nEdges)
edgeFlow=zeros(nEdges,size(M,2));
for ii=1:numel(classIDs)
    j=classIDs(ii); paths=classPaths{j}; w=classPathWeights{j};
    for k=1:numel(paths)
        e=unique(paths{k}(:)); e=e(e>=1&e<=nEdges); if isempty(e), continue; end
        edgeFlow(e,:)=edgeFlow(e,:)+w(k)*M(j,:);
    end
end
end
function geom=directedEdgeGeometryLocal(G,edgeIDs)
[s,t]=edgeEndpointIndicesLocal(G,edgeIDs); lat=double(G.Nodes.Lat); lon=double(G.Nodes.Lon);
geom.s=s(:); geom.t=t(:); geom.lat1=lat(s); geom.lon1=lon(s); geom.lat2=lat(t); geom.lon2=lon(t); geom.lengthKm=haversineVectorKmLocal(geom.lat1,geom.lon1,geom.lat2,geom.lon2);
end
function out=subsetEdgeGeometryLocal(g,mask)
f=fieldnames(g); for k=1:numel(f), out.(f{k})=g.(f{k})(mask); end
end
function g=emptyDisplayGeomLocal()
g=struct();
g.lat1=zeros(0,1); g.lon1=zeros(0,1); g.lat2=zeros(0,1); g.lon2=zeros(0,1);
g.pathLengthKm=zeros(0,1); g.originalEdgeCount=zeros(0,1);
g.latPath=cell(0,1); g.lonPath=cell(0,1); g.cumKm=cell(0,1);
end
function g=displayGeomFromOriginalLocal(edgeGeom,idx)
idx=idx(:); g=emptyDisplayGeomLocal(); if isempty(idx), return; end
n=numel(idx);
g.lat1=edgeGeom.lat1(idx); g.lon1=edgeGeom.lon1(idx); g.lat2=edgeGeom.lat2(idx); g.lon2=edgeGeom.lon2(idx);
g.pathLengthKm=edgeGeom.lengthKm(idx); g.originalEdgeCount=ones(n,1);
g.latPath=cell(n,1); g.lonPath=cell(n,1); g.cumKm=cell(n,1);
for q=1:n
    k=idx(q); L=max(edgeGeom.lengthKm(k),eps);
    g.latPath{q}=[edgeGeom.lat1(k);edgeGeom.lat2(k)];
    g.lonPath{q}=[edgeGeom.lon1(k);edgeGeom.lon2(k)];
    g.cumKm{q}=[0;L];
end
end
function [latLine,lonLine]=displaySegmentsLatLonLocal(g,idx)
% One NaN-separated polyline per merged display segment.  Intermediate road
% nodes are preserved, so long merged chains are still drawn along the road.
idx=idx(:);
if isempty(idx), latLine=zeros(0,1); lonLine=zeros(0,1); return; end
if isfield(g,'latPath') && numel(g.latPath)>=max(idx)
    nTot=sum(cellfun(@numel,g.latPath(idx)))+numel(idx);
    latLine=nan(nTot,1); lonLine=nan(nTot,1); k=0;
    for ii=1:numel(idx)
        q=idx(ii); la=g.latPath{q}(:); lo=g.lonPath{q}(:); m=numel(la);
        latLine(k+(1:m))=la; lonLine(k+(1:m))=lo; k=k+m+1;
    end
else
    n=numel(idx); latLine=nan(3*n,1); lonLine=nan(3*n,1);
    latLine(1:3:end)=g.lat1(idx); latLine(2:3:end)=g.lat2(idx);
    lonLine(1:3:end)=g.lon1(idx); lonLine(2:3:end)=g.lon2(idx);
end
end
function [g,segFlow,carWeight,info]=mergeDisplayEdgesLocal(edgeGeom,rate,flowTol,angleTolDeg,maxChainEdges)
% Merge only unambiguous directed chains.  A chain may continue through a
% node only when the positive-flow subgraph has exactly one incoming and one
% outgoing edge there, adjacent flows are similar, and the turn is small.
rate=double(rate(:)); pos=find(isfinite(rate)&rate>0); g=emptyDisplayGeomLocal(); segFlow=zeros(0,1); carWeight=zeros(0,1);
info=struct('originalPositiveEdges',numel(pos),'displaySegments',0,'mergedEdgesSaved',0);
if isempty(pos), return; end
s=edgeGeom.s(:); t=edgeGeom.t(:); nAll=numel(rate); isPos=false(nAll,1); isPos(pos)=true;
maxNode=max([s(pos);t(pos)]); inCount=accumarray(t(pos),1,[maxNode 1]); outCount=accumarray(s(pos),1,[maxNode 1]);
% Unique successor / predecessor local edge indices for unambiguous nodes.
succ=zeros(nAll,1); pred=zeros(nAll,1);
outEdge=zeros(maxNode,1); inEdge=zeros(maxNode,1);
for kk=reshape(pos,1,[])
    if outCount(s(kk))==1, outEdge(s(kk))=kk; end
    if inCount(t(kk))==1, inEdge(t(kk))=kk; end
end
for kk=reshape(pos,1,[])
    node=t(kk); if node<=maxNode && inCount(node)==1 && outCount(node)==1, succ(kk)=outEdge(node); end
    node=s(kk); if node<=maxNode && inCount(node)==1 && outCount(node)==1, pred(kk)=inEdge(node); end
end
visited=false(nAll,1); chains={};
% Prefer true chain starts; cycles are handled in the second pass.
starts=zeros(0,1);
for kk=reshape(pos,1,[])
    pp=pred(kk);
    if pp==0 || ~canMergePairLocal(edgeGeom,rate,pp,kk,flowTol,angleTolDeg), starts(end+1,1)=kk; end %#ok<AGROW>
end
order=[starts; setdiff(pos,starts,'stable')];
for ii=1:numel(order)
    first=order(ii); if visited(first), continue; end
    chain=first; visited(first)=true; cur=first;
    while numel(chain)<maxChainEdges
        nx=succ(cur); if nx==0 || visited(nx), break; end
        if ~canMergePairLocal(edgeGeom,rate,cur,nx,flowTol,angleTolDeg), break; end
        % Prevent gradual bends from collapsing into a poor straight chord.
        if ~chainChordCompatibleLocal(edgeGeom,chain,nx,angleTolDeg), break; end
        chain(end+1)=nx; visited(nx)=true; cur=nx; %#ok<AGROW>
    end
    chains{end+1}=chain; %#ok<AGROW>
end
n=numel(chains); g.lat1=zeros(n,1); g.lon1=zeros(n,1); g.lat2=zeros(n,1); g.lon2=zeros(n,1); g.pathLengthKm=zeros(n,1); g.originalEdgeCount=zeros(n,1);
g.latPath=cell(n,1); g.lonPath=cell(n,1); g.cumKm=cell(n,1);
segFlow=zeros(n,1); carWeight=zeros(n,1);
for q=1:n
    c=chains{q}; L=max(edgeGeom.lengthKm(c),eps); totalL=sum(L);
    latv=[edgeGeom.lat1(c(1)); edgeGeom.lat2(c(:))];
    lonv=[edgeGeom.lon1(c(1)); edgeGeom.lon2(c(:))];
    cum=[0; cumsum(L(:))];
    % Remove any accidental duplicate cumulative-distance vertices.
    keep=[true; diff(cum)>max(eps,1e-12)]; latv=latv(keep); lonv=lonv(keep); cum=cum(keep);
    g.lat1(q)=latv(1); g.lon1(q)=lonv(1); g.lat2(q)=latv(end); g.lon2(q)=lonv(end);
    g.pathLengthKm(q)=cum(end); g.originalEdgeCount(q)=numel(c);
    g.latPath{q}=latv; g.lonPath{q}=lonv; g.cumKm{q}=cum;
    segFlow(q)=sum(rate(c).*L)/totalL;       % length-weighted representative flow
    carWeight(q)=sum(rate(c));               % preserves old per-edge car allocation mass
end
info.displaySegments=n; info.mergedEdgesSaved=numel(pos)-n;
end
function tf=canMergePairLocal(edgeGeom,rate,a,b,flowTol,angleTolDeg)
if a<=0 || b<=0 || edgeGeom.t(a)~=edgeGeom.s(b), tf=false; return; end
fa=rate(a); fb=rate(b); rel=abs(fa-fb)/max([fa fb eps]); if rel>flowTol, tf=false; return; end
v1=localDirectionVectorLocal(edgeGeom.lat1(a),edgeGeom.lon1(a),edgeGeom.lat2(a),edgeGeom.lon2(a));
v2=localDirectionVectorLocal(edgeGeom.lat1(b),edgeGeom.lon1(b),edgeGeom.lat2(b),edgeGeom.lon2(b));
tf=directionAngleDegLocal(v1,v2)<=angleTolDeg;
end
function tf=chainChordCompatibleLocal(edgeGeom,chain,nx,angleTolDeg)
first=chain(1); last=chain(end);
vChord=localDirectionVectorLocal(edgeGeom.lat1(first),edgeGeom.lon1(first),edgeGeom.lat2(nx),edgeGeom.lon2(nx));
vFirst=localDirectionVectorLocal(edgeGeom.lat1(first),edgeGeom.lon1(first),edgeGeom.lat2(first),edgeGeom.lon2(first));
vLast=localDirectionVectorLocal(edgeGeom.lat1(last),edgeGeom.lon1(last),edgeGeom.lat2(nx),edgeGeom.lon2(nx));
tf=directionAngleDegLocal(vChord,vFirst)<=1.5*angleTolDeg && directionAngleDegLocal(vChord,vLast)<=1.5*angleTolDeg;
end
function v=localDirectionVectorLocal(lat1,lon1,lat2,lon2)
latm=0.5*(lat1+lat2); x=(lon2-lon1)*cosd(latm); y=lat2-lat1; n=hypot(x,y); if n<=eps, v=[1 0]; else, v=[x/n y/n]; end
end
function a=directionAngleDegLocal(v1,v2)
d=max(-1,min(1,sum(v1.*v2))); a=acosd(d);
end
function [latLUT,lonLUT]=sampleDisplayGeometryLocal(g,nSamples)
% Resample each merged polyline at equally spaced cumulative-distance phases.
% The timer can then move every car with a pair of vectorized table lookups.
nSeg=numel(g.pathLengthKm); nSamples=max(2,round(nSamples));
latLUT=nan(nSeg,nSamples); lonLUT=nan(nSeg,nSamples);
phase=linspace(0,1,nSamples);
for q=1:nSeg
    L=max(g.pathLengthKm(q),eps); target=phase*L;
    if isfield(g,'cumKm') && numel(g.cumKm)>=q && numel(g.cumKm{q})>=2
        d=g.cumKm{q}(:); la=g.latPath{q}(:); lo=g.lonPath{q}(:);
        [d,ia]=unique(d,'stable'); la=la(ia); lo=lo(ia);
        if numel(d)>=2
            latLUT(q,:)=interp1(d,la,target,'linear','extrap');
            lonLUT(q,:)=interp1(d,lo,target,'linear','extrap');
        else
            latLUT(q,:)=la(1); lonLUT(q,:)=lo(1);
        end
    else
        latLUT(q,:)=g.lat1(q)+(g.lat2(q)-g.lat1(q))*phase;
        lonLUT(q,:)=g.lon1(q)+(g.lon2(q)-g.lon1(q))*phase;
    end
end
end

function counts=proportionalIntegerCountsLocal(weights,totalCount,ensureCoverage)
w=max(double(weights(:)),0); counts=zeros(size(w)); if totalCount<=0||sum(w)<=0, return; end
pos=w>0; base=zeros(size(w));
if ensureCoverage, base(pos)=1; totalCount=max(totalCount,nnz(pos)); end
rem=totalCount-sum(base); counts=base; if rem<=0, return; end
q=rem*w/sum(w); a=floor(q); counts=counts+a; left=rem-sum(a); if left>0, frac=q-a; [~,ord]=sort(frac,'descend'); counts(ord(1:left))=counts(ord(1:left))+1; end
end
function [R,resultFile]=loadResultsLocal(input)
resultFile=''; if ischar(input)||isstring(input), resultFile=char(string(input)); S=load(resultFile); if isfield(S,'results'), R=S.results; else, fn=fieldnames(S); if numel(fn)==1&&isstruct(S.(fn{1})), R=S.(fn{1}); else, R=S; end, end; elseif isstruct(input), R=input; else, error('animateParmaTrafficFlows:Input','Input must be results struct or MAT filename.'); end
end
function [latLine,lonLine]=edgeSegmentsLatLonLocal(G,edgeIDs)
[s,t]=edgeEndpointIndicesLocal(G,edgeIDs); lat=double(G.Nodes.Lat); lon=double(G.Nodes.Lon); n=numel(edgeIDs); latLine=nan(3*n,1); lonLine=nan(3*n,1); latLine(1:3:end)=lat(s); latLine(2:3:end)=lat(t); lonLine(1:3:end)=lon(s); lonLine(2:3:end)=lon(t);
end
function [s,t]=edgeEndpointIndicesLocal(G,edgeIDs)
ends=G.Edges.EndNodes(edgeIDs,:); if isnumeric(ends), s=ends(:,1); t=ends(:,2); else, [tf,s]=ismember(string(ends(:,1)),string(G.Nodes.Name)); [tf2,t]=ismember(string(ends(:,2)),string(G.Nodes.Name)); if ~all(tf&tf2), error('animateParmaTrafficFlows:Graph','Could not decode graph edge endpoints.'); end, end
end
function label=stratumLabelLocal(T,q)
parts=strings(0,1); vars=T.Properties.VariableNames;
if ismember('DayType',vars), parts(end+1)=string(T.DayType(q)); end
if ismember('YearLabel',vars), parts(end+1)=string(T.YearLabel(q)); elseif ismember('Year',vars), parts(end+1)=string(T.Year(q)); end
if ismember('StartTime',vars)&&ismember('EndTime',vars), parts(end+1)=string(T.StartTime(q))+"-"+string(T.EndTime(q)); end
if ismember('SeasonLabel',vars), parts(end+1)=string(T.SeasonLabel(q)); elseif ismember('SeasonStart',vars)&&ismember('SeasonEnd',vars), parts(end+1)=string(T.SeasonStart(q))+"-"+string(T.SeasonEnd(q)); end
label=char(strjoin(parts,' | '));
end
function tf=isTextSelectionLocal(x), tf=ischar(x)||isstring(x)||iscellstr(x); end
function tf=isScalarLogicalLikeLocal(x), tf=(islogical(x)||isnumeric(x))&&isscalar(x)&&any(x==[0 1]); end
function idx=resolvePoissonSliceLocal(T,wanted)
if ischar(wanted)||isstring(wanted), s=string(wanted); if strcmpi(s,'auto'), idx=1; return; end; z=str2double(s); if ~isfinite(z), error('animateParmaTrafficFlows:PoissonSlice','PoissonSlice must be auto or numeric.'); end; wanted=z; end
wanted=double(wanted); if ismember('Stratum',T.Properties.VariableNames), q=find(double(T.Stratum)==wanted,1); if ~isempty(q), idx=q; return; end, end
if wanted>=1&&wanted<=height(T)&&wanted==floor(wanted), idx=wanted; else, error('animateParmaTrafficFlows:PoissonSlice','Requested slice %g not found.',wanted); end
end
function labels=poissonSliceLabelsLocal(T)
labels=strings(height(T),1); for q=1:height(T), labels(q)=string(sprintf('%s  |  %s',poissonSliceIdLocal(T,q),stratumLabelLocal(T,q))); end
end
function id=poissonSliceIdLocal(T,q)
if ismember('Stratum',T.Properties.VariableNames), id=char(string(T.Stratum(q))); else, id=sprintf('%d',q); end
end
function d=haversineVectorKmLocal(lat1,lon1,lat2,lon2)
R=6371.0088; p1=deg2rad(lat1); p2=deg2rad(lat2); dp=deg2rad(lat2-lat1); dl=deg2rad(lon2-lon1); a=sin(dp/2).^2+cos(p1).*cos(p2).*sin(dl/2).^2; d=2*R*atan2(sqrt(a),sqrt(max(0,1-a)));
end
