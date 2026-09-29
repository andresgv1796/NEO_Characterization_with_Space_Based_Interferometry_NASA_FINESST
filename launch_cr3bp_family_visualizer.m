function handles = launch_cr3bp_family_visualizer(systemName,poFile,overlays,epochUTC,opts)
%LAUNCH_CR3BP_FAMILY_VISUALIZER Fast unified explorer for precomputed families.
%
% The CR3BP orbit/sphere panels are rendered directly from OVERLAYS. The
% asteroid routine is called with RenderPlot=false, so no hidden science
% figures are created/copied. This substantially reduces startup graphics
% overhead for large 10-period families.
%
% Range mode: two-thumb diameter/spin/family selectors.
% Single mode: one-thumb selectors for exploring one value/parent at a time.
% Propagation is always a single horizon control and reveals prefixes of the
% already-computed normal histories; dynamics are never rerun here.

arguments
    systemName (1,1) string
    poFile (1,1) string
    overlays (1,:) struct
    epochUTC (1,1) datetime

    opts.Observer (1,1) string {mustBeMember(opts.Observer,["geocentric","heliocentric"])} = "heliocentric"
    opts.Population (1,1) string {mustBeMember(opts.Population,["MBA","all_asteroids"])} = "all_asteroids"
    opts.CoordinateFrame (1,1) string {mustBeMember(opts.CoordinateFrame,["equatorial","ecliptic","galactic"])} = "ecliptic"
    opts.ShowMilkyWay (1,1) logical = true
    opts.ShowPlanets (1,1) logical = false
    opts.ShowPlanetLabels (1,1) logical = false
    opts.ShowPrimaryBody = []
    opts.ShowSecondaryBody (1,1) logical = true
    opts.Verbose (1,1) logical = false
    % Classification around f_v=0.5. With the default 0.10:
    %   vertical  : f_v >= 0.60
    %   planar    : f_v <= 0.40
    %   undecided : 0.40 < f_v < 0.60
    opts.ModeAmbiguityHalfWidth (1,1) double {mustBeNonnegative} = 0.10
    opts.FigurePosition (1,4) double = [30 30 1780 1040]
end

if opts.ModeAmbiguityHalfWidth >= 0.5
    error('ModeAmbiguityHalfWidth must be smaller than 0.5.');
end

if isempty(overlays), error('overlays must be nonempty.'); end
keep = arrayfun(@(o) isfield(o,'vectorsEclipticJ2000') && ~isempty(o.vectorsEclipticJ2000) && ...
    isfield(o,'orbitSynodic') && ~isempty(o.orbitSynodic) && ...
    isfield(o,'parentRow') && isfinite(o.parentRow),overlays);
overlays = overlays(keep);
if isempty(overlays), error('No usable overlays remain.'); end

nOverlays = numel(overlays);
parentRows = double([overlays.parentRow]);
[~,familyName,~] = fileparts(char(poFile));

verticalFraction = NaN(1,nOverlays);
modeClass = strings(1,nOverlays);
for k=1:nOverlays
    if isfield(overlays(k),'centerVerticalFraction') && ...
            isfinite(overlays(k).centerVerticalFraction)
        verticalFraction(k)=double(overlays(k).centerVerticalFraction);
    end

    if ~isfinite(verticalFraction(k))
        modeClass(k)="undecided";
    elseif verticalFraction(k) >= 0.5 + opts.ModeAmbiguityHalfWidth
        modeClass(k)="vertical";
    elseif verticalFraction(k) <= 0.5 - opts.ModeAmbiguityHalfWidth
        modeClass(k)="planar";
    else
        modeClass(k)="undecided";
    end
end

parentPeriodDays = NaN(1,nOverlays);
availablePeriods = NaN(1,nOverlays);
for k=1:nOverlays
    if isfield(overlays(k),'parentPeriod_days') && isfinite(overlays(k).parentPeriod_days)
        parentPeriodDays(k)=double(overlays(k).parentPeriod_days);
    end
    if isfield(overlays(k),'tau') && ~isempty(overlays(k).tau) && isfinite(overlays(k).parentPeriod_TU)
        availablePeriods(k)=max(double(overlays(k).tau(:)))/double(overlays(k).parentPeriod_TU);
    end
end
maxPropagationPeriods=min(10,min(availablePeriods,[],'omitnan'));
if ~isfinite(maxPropagationPeriods) || maxPropagationPeriods<=0
    error('Could not determine the precomputed propagation horizon.');
end
if maxPropagationPeriods < 10-1e-8
    warning(['Propagation is limited to %.3g T_primary by the supplied overlays. ' ...
        'Use NPrimaryPeriods=10 for the full visualizer range.'],maxPropagationPeriods);
end
initialP=min(1,maxPropagationPeriods);
propLimits=[min(0.25,maxPropagationPeriods) maxPropagationPeriods];

familyLimits=[min(parentRows) max(parentRows)];
if familyLimits(2)<=familyLimits(1), familyLimits(2)=familyLimits(1)+1; end
lastRows=[min(parentRows) max(parentRows)];
lastP=initialP;
lastAsteroidCount=0;

%% Visible shell immediately, before asteroid propagation.
fig=uifigure('Theme','light','Color','w','Name','CR3BP Family Visualizer','Position',opts.FigurePosition);
main=uigridlayout(fig,[4 2]);
main.RowHeight={54,86,'0.92x','1.62x'};
main.ColumnWidth={'1x','1x'};
main.Padding=[12 8 12 10]; main.RowSpacing=7; main.ColumnSpacing=12; main.BackgroundColor='w';

titleBox=uigridlayout(main,[2 1]); titleBox.Layout.Row=1; titleBox.Layout.Column=[1 2];
titleBox.RowHeight={27,23}; titleBox.Padding=[0 0 0 0]; titleBox.RowSpacing=0; titleBox.BackgroundColor='w';
mainTitle=uilabel(titleBox,'Text',sprintf('%s | %s',systemName,familyName), ...
    'FontWeight','bold','FontSize',14,'HorizontalAlignment','center','FontColor','k','BackgroundColor','w');
stateLine=uilabel(titleBox,'Text','Preparing family graphics ...','FontSize',11, ...
    'HorizontalAlignment','center','FontColor','k','BackgroundColor','w');

ctrlGrid=uigridlayout(main,[1 5]); ctrlGrid.Layout.Row=2; ctrlGrid.Layout.Column=[1 2];
ctrlGrid.ColumnWidth={'1x','1x','1x','1x',150}; ctrlGrid.Padding=[2 1 2 1]; ctrlGrid.ColumnSpacing=12; ctrlGrid.BackgroundColor='w';

axOrbit=uiaxes(main); axOrbit.Layout.Row=3; axOrbit.Layout.Column=1;
axSky=uiaxes(main); axSky.Layout.Row=3; axSky.Layout.Column=2;
axAst=uiaxes(main); axAst.Layout.Row=4; axAst.Layout.Column=[1 2];
force_light_axes(axOrbit); force_light_axes(axSky); force_light_axes(axAst);

%% Direct CR3BP rendering.
sys=overlays(1).system; mu=overlays(1).massRatio;
if ~isfinite(mu), mu=sys.mu; end
showPrimary=opts.ShowPrimaryBody;
if isempty(showPrimary), showPrimary=~strcmpi(systemName,'SunEarth'); end
[pName,sName]=body_names(overlays);

hold(axOrbit,'on'); grid(axOrbit,'on'); axis(axOrbit,'equal'); view(axOrbit,38,25);
if showPrimary
    try, plot_body_rotating(axOrbit,pName,[-mu 0 0],LengthUnit_km=sys.Lstar_km,ScaleFactor=1,Resolution=60); catch, end
end
if opts.ShowSecondaryBody
    try, plot_body_rotating(axOrbit,sName,[1-mu 0 0],LengthUnit_km=sys.Lstar_km,ScaleFactor=1,Resolution=60); catch, end
end
xlabel(axOrbit,'$x\,[\mathrm{LU}]$','Interpreter','latex'); ylabel(axOrbit,'$y\,[\mathrm{LU}]$','Interpreter','latex'); zlabel(axOrbit,'$z\,[\mathrm{LU}]$','Interpreter','latex');
title(axOrbit,'Parent periodic orbits','Interpreter','latex');

hold(axSky,'on'); axis(axSky,'equal'); view(axSky,38,24);
[xs,ys,zs]=sphere(36); surf(axSky,xs,ys,zs,'FaceColor',[.86 .86 .86],'FaceAlpha',.055,'EdgeColor',[.72 .72 .72],'EdgeAlpha',.12,'HandleVisibility','off');
plot_sphere_grid_fast(axSky);
xlabel(axSky,sphere_label(opts.CoordinateFrame,'x'),'Interpreter','latex'); ylabel(axSky,sphere_label(opts.CoordinateFrame,'y'),'Interpreter','latex'); zlabel(axSky,sphere_label(opts.CoordinateFrame,'z'),'Interpreter','latex');
title(axSky,"Linear Floquet normal patterns: "+opts.CoordinateFrame+" frame",'Interpreter','none');
xlim(axSky,[-1.05 1.05]); ylim(axSky,[-1.05 1.05]); zlim(axSky,[-1.05 1.05]);

hOrbit=gobjects(1,nOverlays); hOrbitStart=gobjects(1,nOverlays); hSky=gobjects(1,nOverlays); hSkyAnti=gobjects(1,nOverlays); hSkyStart=gobjects(1,nOverlays);
frameFull=cell(1,nOverlays); periodCoord=cell(1,nOverlays); orbLo=NaN(nOverlays,3); orbHi=NaN(nOverlays,3);
for k=1:nOverlays
    c=overlays(k).color; X=double(overlays(k).orbitSynodic(:,1:3));
    hOrbit(k)=plot3(axOrbit,X(:,1),X(:,2),X(:,3),'-','Color',c,'LineWidth',1.45,'DisplayName',char(string(overlays(k).label)),'Tag',sprintf('Orbit_%d',k));
    hOrbitStart(k)=plot3(axOrbit,X(1,1),X(1,2),X(1,3),'o','MarkerSize',4,'MarkerFaceColor',c,'MarkerEdgeColor','k','HandleVisibility','off');
    good=all(isfinite(X),2); if any(good), orbLo(k,:)=min(X(good,:),[],1); orbHi(k,:)=max(X(good,:),[],1); end
    V=transform_ecliptic_vectors_fast(overlays(k).vectorsEclipticJ2000,opts.CoordinateFrame);
    frameFull{k}=V; periodCoord{k}=double(overlays(k).tau(:))/double(overlays(k).parentPeriod_TU);
    nk=prefix_count(periodCoord{k},initialP); VV=V(1:nk,:);
    hSky(k)=plot3(axSky,VV(:,1),VV(:,2),VV(:,3),'-','Color',c,'LineWidth',1.55,'DisplayName',char(string(overlays(k).label)),'Tag',sprintf('Sky_%d',k));
    if isfield(overlays(k),'showAntipode') && overlays(k).showAntipode
        hSkyAnti(k)=plot3(axSky,-VV(:,1),-VV(:,2),-VV(:,3),'--','Color',c,'LineWidth',1.2,'HandleVisibility','off');
    end
    hSkyStart(k)=plot3(axSky,VV(1,1),VV(1,2),VV(1,3),'o','MarkerSize',4,'MarkerFaceColor',c,'MarkerEdgeColor','k','HandleVisibility','off');
end
autoscale_orbit(axOrbit,orbLo,orbHi,true(1,nOverlays));
drawnow limitrate

%% Asteroid DATA ONLY: no temporary giant figures/copies.
stateLine.Text='Computing asteroid snapshot ...'; drawnow
[~,snapshot] = plot_JPL_asteroid_sky(opts.Observer,epochUTC, ...
    Population=opts.Population,CoordinateFrame=opts.CoordinateFrame, ...
    ShowMilkyWay=false,ShowPlanets=false,ShowPlanetLabels=false, ...
    NormalOverlay=struct([]),RenderPlot=false,Verbose=opts.Verbose);

Dfull=double(snapshot.diameter_km(:)); Sfull=double(snapshot.spinRate_revDay(:));
Xfull=double(snapshot.mollweideX(:)); Yfull=double(snapshot.mollweideY(:));
SizeFull=double(snapshot.markerArea_pt2(:)); Cfull=Sfull;
finiteAst=isfinite(Dfull)&isfinite(Sfull)&isfinite(Xfull)&isfinite(Yfull);
Dlim=finite_limits(min(Dfull(finiteAst)),max(Dfull(finiteAst)));
Slim=finite_limits(min(Sfull(finiteAst)),max(Sfull(finiteAst)));
lastD=Dlim; lastS=Slim; lastAsteroidCount=nnz(finiteAst);

%% Direct Mollweide rendering.
hold(axAst,'on'); axis(axAst,'off');
draw_mollweide_grid_fast(axAst,180);
if opts.ShowMilkyWay, draw_milkyway_fast(axAst,opts.CoordinateFrame,180); end
draw_ecliptic_fast(axAst,opts.CoordinateFrame,180);
hAst=scatter(axAst,Xfull,Yfull,SizeFull,Cfull,'filled','MarkerEdgeColor',[.15 .15 .15],'MarkerEdgeAlpha',.15,'Tag','Asteroids');
colormap(axAst,turbo(256)); cb=colorbar(axAst,'eastoutside'); cb.Label.String='Spin rate [rev/day]'; cb.Color='k';

hMoll=gobjects(1,nOverlays); hMollAnti=gobjects(1,nOverlays); hMollStart=gobjects(1,nOverlays); hMollEnd=gobjects(1,nOverlays); hMollLabel=gobjects(1,nOverlays);
mollFull=cell(1,nOverlays); mollAntiFull=cell(1,nOverlays);
for k=1:nOverlays
    V=frameFull{k}; [lon,lat]=unit_to_lonlat(V); [xx,yy]=mollweide_fast(lon,lat,180); [xx,yy]=break_seams(xx,yy);
    mollFull{k}=struct('X',xx(:),'Y',yy(:));
    Va=-V; [lonA,latA]=unit_to_lonlat(Va); [xa,ya]=mollweide_fast(lonA,latA,180); [xa,ya]=break_seams(xa,ya); mollAntiFull{k}=struct('X',xa(:),'Y',ya(:));
    nk=prefix_count(periodCoord{k},initialP); c=overlays(k).color;
    hMoll(k)=plot(axAst,xx(1:nk),yy(1:nk),'-','Color',c,'LineWidth',1.35,'HandleVisibility','off');
    if isfield(overlays(k),'showAntipode') && overlays(k).showAntipode
        hMollAnti(k)=plot(axAst,xa(1:nk),ya(1:nk),'--','Color',c,'LineWidth',1.0,'HandleVisibility','off');
    end
    fi=find(isfinite(xx(1:nk))&isfinite(yy(1:nk))); if isempty(fi), fi=1; end
    hMollStart(k)=plot(axAst,xx(fi(1)),yy(fi(1)),'o','MarkerSize',4,'MarkerFaceColor',c,'MarkerEdgeColor','k','HandleVisibility','off');
    hMollEnd(k)=plot(axAst,xx(fi(end)),yy(fi(end)),'s','MarkerSize',4,'MarkerFaceColor',c,'MarkerEdgeColor','k','HandleVisibility','off');
    im=fi(max(1,round(numel(fi)/2))); hMollLabel(k)=text(axAst,xx(im)+.03,yy(im)+.02,char(string(overlays(k).label)),'Color',c,'FontSize',8,'Interpreter','none','Visible','off');
end
% Fixed projection box avoids the initial off-center appearance seen before zooming.
axis(axAst,'equal'); xlim(axAst,[-2.92 2.92]); ylim(axAst,[-1.48 1.48]);
axAst.DataAspectRatio=[1 1 1]; axAst.PlotBoxAspectRatio=[2 1 1]; axAst.PositionConstraint='outerposition';
title(axAst,sprintf('%s asteroid sky | %s',opts.Observer,opts.CoordinateFrame),'Interpreter','none');
drawnow
% Reapply after layout/colorbar settles.
xlim(axAst,[-2.92 2.92]); ylim(axAst,[-1.48 1.48]); axAst.PlotBoxAspectRatio=[2 1 1];

%% Controls: range + single versions occupy the same slots.
[dRange,dSingle,lD]=make_dual_slider(ctrlGrid,1,'Asteroid diameter [km]',Dlim,Dlim,median(Dlim),false);
[sRange,sSingle,lS]=make_dual_slider(ctrlGrid,2,'Spin rate / color [rev/day]',Slim,Slim,median(Slim),false);
[fRange,fSingle,lF]=make_dual_slider(ctrlGrid,3,'Parent family row',familyLimits,lastRows,nearest(parentRows,median(lastRows)),true);
[pSlider,lP]=make_value_slider(ctrlGrid,4,'Propagation horizon',propLimits,initialP);

actions=uigridlayout(ctrlGrid,[6 1]); actions.Layout.Column=5; actions.RowHeight={17,25,17,25,17,31}; actions.Padding=[0 0 0 0]; actions.RowSpacing=1; actions.BackgroundColor='w';
uilabel(actions,'Text','Selection mode','HorizontalAlignment','center','FontColor','k','BackgroundColor','w');
mode=uidropdown(actions,'Items',{'Range','Single'},'Value','Range');
uilabel(actions,'Text','Center character','HorizontalAlignment','center','FontColor','k','BackgroundColor','w');
modeFilter=uidropdown(actions,'Items',{'All','Vertical','Planar','Undecided'},'Value','All');
uilabel(actions,'Text','Export current state','HorizontalAlignment','center','FontColor','k','BackgroundColor','w');
saveButton=uibutton(actions,'Text','Save PNG','FontWeight','bold');

lastP=initialP; singleMode=false; savedD=Dlim; savedS=Slim; savedRows=lastRows;
lastDrag=tic; minPeriod=.075;

% Range callbacks.
dRange.ValueChangingFcn=@(~,e) drag_range('D',e.Value); sRange.ValueChangingFcn=@(~,e) drag_range('S',e.Value); fRange.ValueChangingFcn=@(~,e) drag_range('F',e.Value);
dRange.ValueChangedFcn=@(h,~) finish_range('D',h.Value); sRange.ValueChangedFcn=@(h,~) finish_range('S',h.Value); fRange.ValueChangedFcn=@(h,~) finish_range('F',h.Value);
% Single callbacks.
dSingle.ValueChangingFcn=@(~,e) drag_single('D',e.Value); sSingle.ValueChangingFcn=@(~,e) drag_single('S',e.Value); fSingle.ValueChangingFcn=@(~,e) drag_single('F',e.Value);
dSingle.ValueChangedFcn=@(h,~) finish_single('D',h.Value); sSingle.ValueChangedFcn=@(h,~) finish_single('S',h.Value); fSingle.ValueChangedFcn=@(h,~) finish_single('F',h.Value);
pSlider.ValueChangingFcn=@(~,e) drag_prop(e.Value); pSlider.ValueChangedFcn=@(h,~) finish_prop(h.Value);
mode.ValueChangedFcn=@(~,~) switch_mode();
modeFilter.ValueChangedFcn=@(~,~) center_filter_changed();
saveButton.ButtonPushedFcn=@(~,~) save_png();

update_all(true);

handles=struct('figure',fig,'axesOrbit',axOrbit,'axesSky',axSky,'axesMollweide',axAst, ...
    'snapshot',snapshot,'overlays',overlays,'mode',mode,'modeFilter',modeFilter, ...
    'saveButton',saveButton,'update',@()update_all(true));

    function center_filter_changed()
        update_family(true);
        update_prop();
        update_labels();
        update_state();
        drawnow limitrate
    end

    function switch_mode()
        singleMode=strcmp(mode.Value,'Single');
        if singleMode
            savedD=lastD; savedS=lastS; savedRows=lastRows;
            dv=nearest(Dfull(finiteAst),mean(lastD)); sv=nearest(Sfull(finiteAst),mean(lastS)); fv=nearest(parentRows,mean(lastRows));
            dSingle.Value=dv; sSingle.Value=sv; fSingle.Value=fv;
            set_visible(dRange,false); set_visible(sRange,false); set_visible(fRange,false); set_visible(dSingle,true); set_visible(sSingle,true); set_visible(fSingle,true);
            lastD=[dv dv]; lastS=[sv sv]; lastRows=[fv fv];
        else
            set_visible(dSingle,false); set_visible(sSingle,false); set_visible(fSingle,false); set_visible(dRange,true); set_visible(sRange,true); set_visible(fRange,true);
            dRange.Value=savedD; sRange.Value=savedS; fRange.Value=savedRows; lastD=savedD; lastS=savedS; lastRows=savedRows;
        end
        update_all(true);
    end

    function drag_range(kind,val)
        if toc(lastDrag)<minPeriod, return, end; lastDrag=tic;
        if kind=='D', lastD=sort(double(val)); elseif kind=='S', lastS=sort(double(val)); else, lastRows=snap_row_range(val,parentRows); end
        update_kind(kind,false); drawnow limitrate nocallbacks
    end
    function finish_range(kind,val)
        if kind=='D', lastD=sort(double(val)); savedD=lastD; elseif kind=='S', lastS=sort(double(val)); savedS=lastS; else, lastRows=snap_row_range(val,parentRows); fRange.Value=lastRows; savedRows=lastRows; end
        update_kind(kind,true); drawnow limitrate
    end
    function drag_single(kind,val)
        if toc(lastDrag)<minPeriod, return, end; lastDrag=tic;
        if kind=='D', v=nearest(Dfull(finiteAst),val); lastD=[v v]; elseif kind=='S', v=nearest(Sfull(finiteAst),val); lastS=[v v]; else, v=nearest(parentRows,val); lastRows=[v v]; end
        update_kind(kind,false); drawnow limitrate nocallbacks
    end
    function finish_single(kind,val)
        if kind=='D', v=nearest(Dfull(finiteAst),val); dSingle.Value=v; lastD=[v v]; elseif kind=='S', v=nearest(Sfull(finiteAst),val); sSingle.Value=v; lastS=[v v]; else, v=nearest(parentRows,val); fSingle.Value=v; lastRows=[v v]; end
        update_kind(kind,true); drawnow limitrate
    end
    function drag_prop(val)
        lP.Text=prop_text(val,physical_time(val,lastRows,parentRows,parentPeriodDays));
        if toc(lastDrag)<minPeriod, return, end; lastDrag=tic; lastP=clamp(val,propLimits); update_prop(); drawnow limitrate nocallbacks
    end
    function finish_prop(val), lastP=clamp(val,propLimits); pSlider.Value=lastP; update_prop(); drawnow limitrate; end
    function update_kind(kind,relegend)
        if kind=='D' || kind=='S', update_ast(); else, update_family(relegend); update_prop(); end
        update_labels(); update_state();
    end
    function update_all(relegend), update_ast(); update_family(relegend); update_prop(); update_labels(); update_state(); end
    function update_ast()
        if singleMode
            tolD=max(1e-12,32*eps(max(1,abs(lastD(1))))); tolS=max(1e-12,32*eps(max(1,abs(lastS(1)))));
            mask=finiteAst & abs(Dfull-lastD(1))<=tolD & abs(Sfull-lastS(1))<=tolS;
        else
            mask=finiteAst&Dfull>=lastD(1)&Dfull<=lastD(2)&Sfull>=lastS(1)&Sfull<=lastS(2);
        end
        hAst.XData=Xfull(mask).'; hAst.YData=Yfull(mask).'; hAst.SizeData=SizeFull(mask).'; hAst.CData=Cfull(mask).'; lastAsteroidCount=nnz(mask);
        lo=lastS(1); hi=lastS(2); if hi<=lo, hi=lo+max(1e-6,abs(lo)*1e-6); end; clim(axAst,[lo hi]);
    end
    function selected=selection_mask()
        selected=parentRows>=lastRows(1)&parentRows<=lastRows(2);
        switch string(modeFilter.Value)
            case "Vertical"
                selected=selected & modeClass=="vertical";
            case "Planar"
                selected=selected & modeClass=="planar";
            case "Undecided"
                selected=selected & modeClass=="undecided";
        end
    end
    function update_family(relegend)
        sel=selection_mask(); ns=nnz(sel);
        for j=1:nOverlays
            vv=onoff(sel(j)); hOrbit(j).Visible=vv; hOrbitStart(j).Visible=vv; hSky(j).Visible=vv; hSkyStart(j).Visible=vv;
            if isgraphics(hSkyAnti(j)), hSkyAnti(j).Visible=vv; end
            hMoll(j).Visible=vv; hMollStart(j).Visible=vv; hMollEnd(j).Visible=vv;
            if isgraphics(hMollAnti(j)), hMollAnti(j).Visible=vv; end
            hMollLabel(j).Visible=onoff(sel(j)&&ns<=5);
        end
        autoscale_orbit(axOrbit,orbLo,orbHi,sel);
        if relegend, update_leg(axOrbit,hOrbit,sel); update_leg(axSky,hSky,sel); end
    end
    function update_prop()
        sel=find(selection_mask());
        for j=sel
            nk=prefix_count(periodCoord{j},lastP); V=frameFull{j};
            hSky(j).XData=V(1:nk,1); hSky(j).YData=V(1:nk,2); hSky(j).ZData=V(1:nk,3);
            if isgraphics(hSkyAnti(j)), hSkyAnti(j).XData=-V(1:nk,1); hSkyAnti(j).YData=-V(1:nk,2); hSkyAnti(j).ZData=-V(1:nk,3); end
            M=mollFull{j}; hMoll(j).XData=M.X(1:nk); hMoll(j).YData=M.Y(1:nk);
            if isgraphics(hMollAnti(j)), A=mollAntiFull{j}; hMollAnti(j).XData=A.X(1:nk); hMollAnti(j).YData=A.Y(1:nk); end
            fi=find(isfinite(M.X(1:nk))&isfinite(M.Y(1:nk))); if ~isempty(fi), hMollEnd(j).XData=M.X(fi(end)); hMollEnd(j).YData=M.Y(fi(end)); im=fi(max(1,round(numel(fi)/2))); hMollLabel(j).Position(1:2)=[M.X(im)+.03 M.Y(im)+.02]; end
        end
    end
    function update_labels()
        if singleMode, lD.Text=sprintf('%.3g',lastD(1)); lS.Text=sprintf('%.3g',lastS(1)); lF.Text=sprintf('%d',round(lastRows(1))); else, lD.Text=range_text(lastD); lS.Text=range_text(lastS); lF.Text=sprintf('[%d, %d]',round(lastRows)); end
        lP.Text=prop_text(lastP,physical_time(lastP,lastRows,parentRows,parentPeriodDays));
    end
    function update_state()
        sel=selection_mask(); familyText=compact_family(sel,parentRows); phys=physical_time(lastP,lastRows,parentRows,parentPeriodDays);
        if singleMode, ds=sprintf('D=%.3g km | spin=%.3g rev/day',lastD(1),lastS(1)); else, ds=sprintf('D=%.3g--%.3g km | spin=%.3g--%.3g rev/day',lastD(1),lastD(2),lastS(1),lastS(2)); end
        if string(modeFilter.Value)=="All"
            modeText=sprintf('centers V/P/U=%d/%d/%d', ...
                nnz(modeClass=="vertical"),nnz(modeClass=="planar"),nnz(modeClass=="undecided"));
        else
            modeText=lower(string(modeFilter.Value))+" centers";
        end
        stateLine.Text=sprintf('%s | %s | %s | horizon %.2f Tp (%s) | %d asteroids', ...
            ds,familyText,modeText,lastP,physical_text(phys),lastAsteroidCount);
        fig.Name=sprintf('CR3BP Family Visualizer | %s | %.2f Tp',familyName,lastP);
    end
    function save_png()
        update_all(true); drawnow; sel=selection_mask(); phys=physical_time(lastP,lastRows,parentRows,parentPeriodDays);
        default=sprintf('%s_%s_%s_Tp%.2f_%s.png',systemName,familyName,filename_family(sel,parentRows),lastP,physical_text(phys)); default=regexprep(default,'[^A-Za-z0-9_.-]','_');
        [f,pth]=uiputfile('*.png','Save current visualization',default); if isequal(f,0), return, end
        try, exportapp(fig,fullfile(pth,f)); catch, exportgraphics(fig,fullfile(pth,f),'Resolution',300); end
    end
end

%% UI helpers
function [r,s,l]=make_dual_slider(parent,col,titleText,lims,rval,sval,isInt)
box=uigridlayout(parent,[2 2]); box.Layout.Column=col; box.RowHeight={24,'1x'}; box.ColumnWidth={'1x',125}; box.Padding=[0 0 0 0]; box.RowSpacing=2; box.ColumnSpacing=4; box.BackgroundColor='w';
lab=uilabel(box,'Text',titleText,'FontWeight','bold','FontColor','k','BackgroundColor','w'); lab.Layout.Row=1; lab.Layout.Column=1;
l=uilabel(box,'Text',range_text(rval),'HorizontalAlignment','right','FontColor','k','BackgroundColor','w'); l.Layout.Row=1; l.Layout.Column=2;
r=uislider(box,"range",'Limits',double(lims),'Value',double(rval),'MinorTicks',[]); r.Layout.Row=2; r.Layout.Column=[1 2];
s=uislider(box,'Limits',double(lims),'Value',double(sval),'MinorTicks',[],'Visible','off'); s.Layout.Row=2; s.Layout.Column=[1 2];
if isInt
    ticks=unique(round(linspace(lims(1),lims(2),min(7,max(2,round(lims(2)-lims(1)+1)))))); r.MajorTicks=ticks; r.Step=1; s.MajorTicks=ticks;
else, r.MajorTicks=[]; s.MajorTicks=[]; end
end
function [s,l]=make_value_slider(parent,col,titleText,lims,val)
box=uigridlayout(parent,[2 2]); box.Layout.Column=col; box.RowHeight={24,'1x'}; box.ColumnWidth={'1x',155}; box.Padding=[0 0 0 0]; box.RowSpacing=2; box.ColumnSpacing=4; box.BackgroundColor='w';
lab=uilabel(box,'Text',titleText,'FontWeight','bold','FontColor','k','BackgroundColor','w'); lab.Layout.Row=1; lab.Layout.Column=1;
l=uilabel(box,'Text',sprintf('%.2f Tp',val),'HorizontalAlignment','right','FontColor','k','BackgroundColor','w'); l.Layout.Row=1; l.Layout.Column=2;
s=uislider(box,'Limits',double(lims),'Value',double(val),'MinorTicks',[]); s.Layout.Row=2; s.Layout.Column=[1 2]; s.MajorTicks=unique([lims(1) 1:floor(lims(2)) lims(2)]);
end

%% Plot helpers
function draw_mollweide_grid_fast(ax,center)
t=linspace(0,2*pi,500); plot(ax,2*sqrt(2)*cos(t),sqrt(2)*sin(t),'-','Color',[.25 .25 .25],'LineWidth',.8,'HandleVisibility','off');
for la=-60:30:60, lo=linspace(0,360,721); [x,y]=mollweide_fast(lo,la+zeros(size(lo)),center); plot(ax,x,y,':','Color',[.72 .72 .72],'LineWidth',.5,'HandleVisibility','off'); end
for lo0=0:30:330, la=linspace(-89.5,89.5,361); [x,y]=mollweide_fast(lo0+zeros(size(la)),la,center); plot(ax,x,y,':','Color',[.78 .78 .78],'LineWidth',.45,'HandleVisibility','off'); end
end
function draw_milkyway_fast(ax,frame,center)
try, [~,n,s]=milkyway(); catch, return, end
for B={n,s}
    q=B{1}; V=lonlat_to_unit(q(:,1),q(:,2)); V=transform_equatorial_vectors_fast(V,frame); [lo,la]=unit_to_lonlat(V); [x,y]=mollweide_fast(lo,la,center); [x,y]=break_seams(x,y); plot(ax,x,y,'-','Color',[.45 .45 .45],'LineWidth',.7,'HandleVisibility','off');
end
end
function draw_ecliptic_fast(ax,frame,center)
lo=linspace(0,360,721); V=lonlat_to_unit(lo,zeros(size(lo))); V=transform_ecliptic_vectors_fast(V,frame); [L,B]=unit_to_lonlat(V); [x,y]=mollweide_fast(L,B,center); [x,y]=break_seams(x,y); plot(ax,x,y,'-','Color',[.35 .35 .35],'LineWidth',.8,'HandleVisibility','off');
end
function [x,y]=mollweide_fast(lon,lat,center)
lr=mod(double(lon)-center+180,360)-180; lam=deg2rad(lr); phi=deg2rad(double(lat)); th=phi;
for ii=1:8, den=2+2*cos(2*th); bad=abs(den)<1e-12; step=zeros(size(th)); step(~bad)=(2*th(~bad)+sin(2*th(~bad))-pi*sin(phi(~bad)))./den(~bad); th=th-step; end
x=2*sqrt(2)/pi*lam.*cos(th); y=sqrt(2)*sin(th);
end
function [x,y]=break_seams(x,y), x=x(:); y=y(:); j=[false;abs(diff(x))>1.6]; x(j)=NaN; y(j)=NaN; end
function V=lonlat_to_unit(lon,lat), lon=deg2rad(double(lon(:))); lat=deg2rad(double(lat(:))); V=[cos(lat).*cos(lon) cos(lat).*sin(lon) sin(lat)]; end
function [lon,lat]=unit_to_lonlat(V), V=normalize_rows(V); lon=mod(rad2deg(atan2(V(:,2),V(:,1))),360); lat=rad2deg(asin(V(:,3))); end
function W=transform_ecliptic_vectors_fast(V,frame)
V=normalize_rows(double(V)); if frame=="ecliptic", W=V; elseif frame=="equatorial", W=ecl2eq(V); else, W=eq2gal(ecl2eq(V)); end; W=normalize_rows(W);
end
function W=transform_equatorial_vectors_fast(V,frame)
V=normalize_rows(double(V)); if frame=="equatorial", W=V; elseif frame=="ecliptic", e=deg2rad(84381.448/3600); R=[1 0 0;0 cos(e) sin(e);0 -sin(e) cos(e)]; W=(R*V.').'; else, W=eq2gal(V); end; W=normalize_rows(W);
end
function W=ecl2eq(V), e=deg2rad(84381.448/3600); R=[1 0 0;0 cos(e) -sin(e);0 sin(e) cos(e)]; W=(R*V.').'; end
function W=eq2gal(V), R=[-.0548755604162154 -.873437090234885 -.4838350155487132;.4941094278755837 -.4448296299600112 .7469822444972189;-.8676661490190047 -.1980763734312015 .4559837761750669]; W=(R*V.').'; end
function V=normalize_rows(V), n=vecnorm(V,2,2); V=V./n; end
function plot_sphere_grid_fast(ax)
p=linspace(0,2*pi,241); plot3(ax,cos(p),sin(p),zeros(size(p)),':','Color',[.5 .5 .5],'HandleVisibility','off'); la=linspace(-pi/2,pi/2,121);
for d=0:45:315, l=deg2rad(d); plot3(ax,cos(la)*cos(l),cos(la)*sin(l),sin(la),':','Color',[.78 .78 .78],'LineWidth',.45,'HandleVisibility','off'); end
end
function s=sphere_label(frame,a), if frame=="ecliptic", if a=='x',s='$\hat e_{\lambda,0}$';elseif a=='y',s='$\hat e_{\lambda,90}$';else,s='$\hat e_{\beta,+}$';end; else, s=['$\hat ' a '$']; end; end

%% State/helpers
function n=prefix_count(p,T), n=find(p<=T+64*eps(max(1,T)),1,'last'); if isempty(n),n=1;end; end
function [p,s]=body_names(o), p="Primary";s="Secondary"; if isfield(o(1),'embedding'), E=o(1).embedding; if isfield(E,'primaryName'),p=string(E.primaryName);end;if isfield(E,'secondaryName'),s=string(E.secondaryName);end;end; end
function force_light_axes(a), a.Color='w';a.XColor='k';a.YColor='k';a.ZColor='k'; end
function autoscale_orbit(ax,lo,hi,sel), idx=find(sel); if isempty(idx),return,end; L=min(lo(idx,:),[],1,'omitnan');H=max(hi(idx,:),[],1,'omitnan'); if any(~isfinite(L))||any(~isfinite(H)),return,end; sp=H-L; ref=max(sp); if ref<=0,ref=1e-3;end; for j=1:3, pad=.08*(sp(j)+(sp(j)<=1e-12)*ref); L(j)=L(j)-pad;H(j)=H(j)+pad;end;xlim(ax,[L(1) H(1)]);ylim(ax,[L(2) H(2)]);zlim(ax,[L(3) H(3)]); end
function update_leg(ax,h,sel), try,legend(ax,'off');catch,end; idx=find(sel); if numel(idx)<=5&&~isempty(idx), legend(ax,h(idx),'Location','best','Interpreter','none');end; end
function r=snap_row_range(v,rows), v=sort(double(v)); r=[nearest(rows,v(1)) nearest(rows,v(2))];r=sort(r); end
function v=nearest(a,x), a=double(a(:)); [~,i]=min(abs(a-double(x)));v=a(i); end
function lim=finite_limits(a,b), if b<=a,b=a+max(1,abs(a)*.01);end;lim=[a b];end
function x=clamp(x,l),x=min(max(double(x),l(1)),l(2));end
function t=range_text(v),v=double(v);if abs(v(2)-v(1))<1e-12,t=sprintf('%.3g',v(1));else,t=sprintf('[%.3g, %.3g]',v(1),v(2));end;end
function t=compact_family(sel,rows), rr=rows(sel);n=numel(rr);if n==0,t='0 parents';elseif n==1,t=sprintf('row %d',rr);else,t=sprintf('%d parents, rows %d--%d',n,min(rr),max(rr));end;end
function t=filename_family(sel,rows), rr=rows(sel);if isempty(rr),t='parents0';elseif numel(rr)==1,t=sprintf('row%d',rr);else,t=sprintf('rows%d-%d_N%d',min(rr),max(rr),numel(rr));end;end
function d=physical_time(P,rowRange,rows,periodDays), sel=rows>=rowRange(1)&rows<=rowRange(2); x=P*periodDays(sel);x=x(isfinite(x));if isempty(x),d=[NaN NaN];else,d=[min(x) max(x)];end;end
function t=physical_text(d), if any(~isfinite(d)),t='time unavailable';return,end;if max(d)>=365.25,v=d/365.25;u='yr';else,v=d;u='days';end;if abs(diff(v))<=1e-8*max(1,max(v)),t=sprintf('%.3g %s',v(1),u);else,t=sprintf('%.3g--%.3g %s',v(1),v(2),u);end;end
function t=prop_text(P,d),t=sprintf('%.2f Tp | %s',P,physical_text(d));end
function set_visible(h,tf),h.Visible=onoff(tf);end
function s=onoff(tf),if tf,s='on';else,s='off';end;end
