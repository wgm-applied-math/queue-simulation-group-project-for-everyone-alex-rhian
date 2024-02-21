% Script that runs a ServiceQueue simulation many times and plots a
% histogram

%% Set up

% Set up to run 100 samples of the queue.
n_samples = 100;

% Each sample is run up to a maximum time of 1000.
max_time = 480;

% Record how many customers are in the system at the end of each sample.
NInSystemSamples = cell([1, n_samples]);

%% Run the queue simulation

% The statistics seem to come out a little weird if the log interval is too
% short, apparently because the log entries are not independent enough.  So
% the log interval should be long enough for several arrival and departure
% events happen.
for sample_num = 1:n_samples
    q = ServiceQueue(LogInterval=10);
    q.schedule_event(Arrival(1, Customer(1)));
    run_until(q, max_time);
    % Pull out samples of the number of customers in the queue system. Each
    % sample run of the queue results in a column of samples of customer
    % counts, because tables like q.Log allow easy extraction of whole
    % columns like this.
    NInSystemSamples{sample_num} = q.Log.NWaiting + q.Log.NInService;
end

% Join all the samples. "vertcat" is short for "vertical concatenate",
% meaning it joins a bunch of arrays vertically, which in this case results
% in one tall column.
NInSystem = vertcat(NInSystemSamples{:});

% MATLAB-ism: When you pull multiple items from a cell array, the result is
% a "comma-separated list" rather than some kind of array.  Thus, the above
% means
%
%    NInSystem = horzcat(NInSystemSamples{1}, NInSystemSamples{2}, ...)
%
% which horizontally concatenates all the lists of numbers in
% NInSystemSamples.
%
% This is roughly equivalent to "splatting" in Python, which looks like
% f(*args).

%% Make a picture
qlength = length(q.Served);

W = zeros(1, qlength);

for n = 1:qlength
    W(1, n) = q.Served{1, n}.DepartureTime - q.Served{1, n}.ArrivalTime;
end

totaltimeinsystem = sum(W)/qlength;


WQ = zeros(1, qlength);

for n = 1:qlength
    WQ(1, n) = q.Served{1, n}.BeginServiceTime - q.Served{1, n}.ArrivalTime;
end

totaltimewaiting = sum(WQ)/qlength;


TotalServed = zeros(1, qlength);

for n = 1:qlength
    TotalServed(1, n) = q.Served{1, n}.DepartureTime - q.Served{1, n}.BeginServiceTime;
end

timeserved = sum(TotalServed)/qlength;


SimulationSteadyStates = [totaltimeinsystem, totaltimewaiting, timeserved];
TheoreticalSteadyStates = [5.361, .3024, 5.0586];

figw = figure();
tw = tiledlayout(figw, 1, 1);
axw = nexttile(tw);
hold(axw, "on");

Wh = histogram(axw, W, Normalization = "probability", BinMethod = "auto");
plot(axw, xline(5.361));

figwq = figure();
twq = tiledlayout(figwq, 1, 1);
axwq = nexttile(twq);
hold(axwq, "on");

Wqh = histogram(axwq, WQ, Normalization = "probability", BinMethod = "auto");
plot(axwq, xline(.3024));

figts = figure();
tts = tiledlayout(figts, 1, 1);
axts = nexttile(tts);
hold(axts, "on");

TSh = histogram(axts, TotalServed, Normalization = "probability", BinMethod = "auto");
plot(axts, xline(5.0586));


% Start with a histogram.  The result is an empirical PDF, that is, the
% area of the bar at horizontal index n is proportional to the fraction of
% samples for which there were n customers in the system.
fig2 = figure();
t2 = tiledlayout(1,1);
ax2 = nexttile(t2);
h = histogram(NInSystem, Normalization = "probability", BinMethod = "integers");


% MATLAB-ism: Once you've created a picture, you can use "hold on" to cause
% further plotting function to work with the same picture rather than
% create a new one.
hold(ax2, "on");
theoryprob = [.40315, .40315, .15118, .03779, .004716];
xvaluesshifted = [0, 1, 2, 3, 4];
plot(ax2, xvaluesshifted, theoryprob, 'o');


%Wt = histogram(totaltimeinsystem, Normalization= "percentage", BinMethod="integers");
%WQt = histogram(totaltimewaiting, Normalization= "percentage", BinMethod="integers");
%WmWQt = histogram(timeserved, Normalization= "percentage", BinMethod="integers");


% For comparison, plot the theoretical results for a M/M/1 queue.
% The agreement isn't all that good unless you run for a long time, say
% max_time = 10,000 units, and LogInterval is large, say 10.
rho = q.ArrivalRate / q.DepartureRate;
P0 = 1 - rho;
nMax = 10;
ns = 0:nMax;
P = zeros([1, nMax+1]);
P(1) = P0;
for n = 1:nMax
    P(1+n) = P0 * rho^n;
end
plot(ax2, ns, P, 'o', MarkerEdgeColor='k', MarkerFaceColor='r');

% This sets some paper-related properties of the figure so that you can
% save it as a PDF and it doesn't fill a whole page.
% gcf is "get current figure handle"
% See https://stackoverflow.com/a/18868933/2407278
fig = gcf;
fig.Units = 'inches';
screenposition = fig.Position;
fig.PaperPosition = [0 0 screenposition(3:4)];
fig.PaperSize = [screenposition(3:4)];