% Inputs:
%   sacFiles    SAC filename, wildcard string, or cell array of filenames.
%
% Name-value parameters:
%   'PeriodBand'    Period band [Tmin Tmax] in seconds, default [1 30].
%                   This means frequency band [1/Tmax, 1/Tmin] Hz.
%   'FreqBand'      Direct frequency band [fmin fmax] in Hz. If provided,
%                   it overrides PeriodBand.
%   'FilterOrder'   Butterworth filter order, default 4.
%   'Normalize'     Normalize each trace by max(abs(trace)), default true.
%   'TimeXLim'      Time axis limit, default [].
%   'FreqXLim'      Frequency axis limit, default [].
%   'Colors'        Cell array or numeric color matrix, default auto.
%   'LineWidth'     Plot line width, default 1.0.
%   'TitleText'     Figure title text, default generated from SAC metadata.
%   'ApplyFilter'   Whether apply bandpass filter, default true.
%
% Output:
%   out.sacData     Original SAC structure array returned by readsac.
%   out.meta        Important SAC metadata.
%   out.waveform    Time, raw data, filtered data.
%   out.spectrum    Frequency and normalized amplitude spectrum.
function out=Fun_PlotWaves(sacFiles, varargin)
p=inputParser;
addRequired(p, 'sacFiles', @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'PeriodBand', [1 30], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'FreqBand', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
addParameter(p, 'FilterOrder', 4, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Normalize', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'TimeXLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
addParameter(p, 'FreqXLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
addParameter(p, 'LineWidth', 1.0, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'TitleText', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'ApplyFilter', true, @(x) islogical(x) || isnumeric(x));
parse(p, sacFiles, varargin{:});

periodBand = sort(p.Results.PeriodBand);
freqBand = p.Results.FreqBand;
filterOrder = p.Results.FilterOrder;
normalizeTrace = logical(p.Results.Normalize);
timeXLim = p.Results.TimeXLim;
freqXLim = p.Results.FreqXLim;
lineWidth = p.Results.LineWidth;
titleText = string(p.Results.TitleText);
applyFilterFlag = logical(p.Results.ApplyFilter);

fileList = Fun_ExpandSacFileList(sacFiles);
if isempty(fileList)
    error('No SAC files found.');
end

sacData = readsac(fileList);
if isempty(sacData)
    error('readsac returned empty data.');
end

if numel(sacData) ~= numel(fileList)
    warning('Number of readable SAC files is %d, but input files are %d.', ...
        numel(sacData), numel(fileList));
end
nSac = numel(sacData);

out = struct();
out.sacData = sacData;
out.meta = struct([]);
out.waveform = struct([]);
out.spectrum = struct([]);

% Figure setting
figure;
tiledlayout(2,1,"TileSpacing","compact","Padding","compact");
nexttile();
ax1=gca;
Fun_SetAxis(ax1);
hold(ax1, 'on');
nexttile();
ax2=gca;
Fun_SetAxis(ax2);
hold(ax2, 'on');
legendNames = strings(nSac, 1);
% Loop the colors
cols=[0.5020    0.5020    0.5020;
      0.0745    0.6235    1;
      0.0314    0.7804    0.4314;
      0.9098    0.0627    0.4039;
      0.4941    0.1843    0.5569;
      1.0000    0.4118    0.1608;
      0     0     0;
      0.6353    0.0784    0.1843];
cols=repmat(cols, 10, 1);
for i = 1:nSac
    s = sacData(i);

    dt = s.DELTA;
    npts = s.NPTS;
    fs = 1 / dt;

    if isnan(dt) || dt <= 0
        error('Invalid DELTA in file: %s', s.FILENAME);
    end

    if isnan(npts) || npts <= 1
        error('Invalid NPTS in file: %s', s.FILENAME);
    end

    if isfield(s, 'B') && ~isnan(s.B)
        bTime = s.B;
    else
        bTime = 0;
    end

    t = bTime + (0:npts-1)' * dt;
    xRaw = double(s.DATA1(:));

    xRaw = detrend(xRaw, 'constant');
    xRaw = detrend(xRaw, 'linear');

    if isempty(freqBand)
        fLow = 1 / periodBand(2);
        fHigh = 1 / periodBand(1);
        bandText = sprintf('bp %.3g-%.3g s', periodBand(1), periodBand(2));
    else
        fb = sort(freqBand);
        fLow = fb(1);
        fHigh = fb(2);
        bandText = sprintf('bp %.3g-%.3g Hz', fLow, fHigh);
    end

    if applyFilterFlag
        xFilt = Fun_ButterPass(xRaw, fs, fLow, fHigh, filterOrder);
    else
        xFilt = xRaw;
        bandText = 'unfiltered';
    end

    if normalizeTrace
        scale = max(abs(xFilt), [], 'omitnan');
        if isempty(scale) || scale == 0 || isnan(scale)
            scale = 1;
        end
        xPlot = xFilt / scale + i;
    else
        xPlot = xFilt + i;
    end

    %     col=getColorLocal(colors,i);
    col=cols(i,:);
    plot(ax1, t, xPlot, ...
        'Color', col, 'LineWidth', lineWidth);
    set(ax1, 'YDir', 'reverse');

    [freq, ampNorm, ampRaw] = Fun_CalFFT_Single(xFilt, fs);
    semilogy(ax2, freq, ampNorm, ...
        'Color', col, 'LineWidth', lineWidth);
    legendNames(i) = getLegendNameLocal(s, i);
    % Set the log-space or not
    ax2.XScale='log';
    ax2.YScale='log';
    ylim(ax2,[1e-8 1]);
    out.meta(i).filename = getFieldOrDefaultLocal(s, 'FILENAME', '');
    out.meta(i).network = strtrim(getFieldOrDefaultLocal(s, 'KNETWK', ''));
    out.meta(i).station = strtrim(getFieldOrDefaultLocal(s, 'KSTNM', ''));
    out.meta(i).component = strtrim(getFieldOrDefaultLocal(s, 'KCMPNM', ''));
    out.meta(i).delta = getFieldOrDefaultLocal(s, 'DELTA', NaN);
    out.meta(i).fs = fs;
    out.meta(i).npts = getFieldOrDefaultLocal(s, 'NPTS', NaN);
    out.meta(i).b = getFieldOrDefaultLocal(s, 'B', NaN);
    out.meta(i).e = getFieldOrDefaultLocal(s, 'E', NaN);
    out.meta(i).a = getFieldOrDefaultLocal(s, 'A', NaN);
    out.meta(i).o = getFieldOrDefaultLocal(s, 'O', NaN);
    out.meta(i).stla = getFieldOrDefaultLocal(s, 'STLA', NaN);
    out.meta(i).stlo = getFieldOrDefaultLocal(s, 'STLO', NaN);
    out.meta(i).stel = getFieldOrDefaultLocal(s, 'STEL', NaN);
    out.meta(i).evla = getFieldOrDefaultLocal(s, 'EVLA', NaN);
    out.meta(i).evlo = getFieldOrDefaultLocal(s, 'EVLO', NaN);
    out.meta(i).evdp = getFieldOrDefaultLocal(s, 'EVDP', NaN);
    out.meta(i).mag = getFieldOrDefaultLocal(s, 'MAG', NaN);
    out.meta(i).dist = getFieldOrDefaultLocal(s, 'DIST', NaN);
    out.meta(i).az = getFieldOrDefaultLocal(s, 'AZ', NaN);
    out.meta(i).baz = getFieldOrDefaultLocal(s, 'BAZ', NaN);
    out.meta(i).gcarc = getFieldOrDefaultLocal(s, 'GCARC', NaN);
    out.meta(i).idep = getFieldOrDefaultLocal(s, 'IDEP', '');
    out.meta(i).iftype = getFieldOrDefaultLocal(s, 'IFTYPE', '');

    out.waveform(i).filename = out.meta(i).filename;
    out.waveform(i).time = t;
    out.waveform(i).raw = xRaw;
    out.waveform(i).filtered = xFilt;
    out.waveform(i).plotData = xPlot;

    out.spectrum(i).filename = out.meta(i).filename;
    out.spectrum(i).freq = freq;
    out.spectrum(i).amplitude = ampRaw;
    out.spectrum(i).amplitudeNorm = ampNorm;
end
xlabel(ax1,'Time (s)');
ylabel(ax1,'Trace');
yticks(ax1,1:numel(sacFiles));
if isempty(timeXLim)
    xlim(ax1, getCommonXLimLocal(out.waveform, 'time'));
else
    xlim(ax1, timeXLim);
end
if titleText == ""
    titleText = buildTitleLocal(out.meta, bandText);
end
xlabel(ax2,'Frequency (Hz)');
ylabel(ax2,'Amplitude');

if isempty(freqXLim)
    maxFs = max([out.meta.fs]);
    xlim(ax2, [1e-2 maxFs/3]);
else
    xlim(ax2, freqXLim);
end

set(gcf,'position',[200,50,900,800]);
end

