%% You may want to change these, first two must be consistent with size of data %
numLags = 4000;
nodesToAvrg = 1000;
FILENAMES = ["timeSrData-tau-0.08.txt", "timeSrData-tau-0.8.txt", "timeSrData-tau-4.0.txt", "timeSrData-tau-40.0.txt"];
taus = ["0.08","0.8","4.0","40.0"];
%% Reading takes time %%
data = {};
for i = 1:length(FILENAMES)
    data{i} = importdata( FILENAMES{i} );
end
%% Read metadata %%
timePoints = [];
for i = 1:length(FILENAMES)
    timePoints(i) = size(data{i},1);
end
totNodes = [];
for i = 1:length(FILENAMES)
    totNodes(i) = size(data{i},2);
end
%% Calculate autocorrelation of selected columns/time-series for each cell in struct %%
ACF_DATA_STRUCT = {};
for i=1:length(FILENAMES)
    autoCorrelations = [];
    signal = data{i};
    for j=1:nodesToAvrg
        ACF = autocorr( signal(:,j), 'numlags', numLags);
        autoCorrelations = [autoCorrelations ACF];
    end
    ACF_DATA_STRUCT{i} = autoCorrelations; %each cell contains a matrix with autocorrelations
end

%% Average them over all points %%
avrgACF_DATA_STRUCT = {};
for i=1:length(FILENAMES)
    avrgACF = [];
    autocorr_data = ACF_DATA_STRUCT{i}; 
    for j=1:numLags
        avrgACF(j) = sum( autocorr_data(j,:) ) / nodesToAvrg;
    end
    avrgACF_DATA_STRUCT{i} = avrgACF; 
end
%% Linear Plot %%
figure(1)
plot(avrgACF_DATA_STRUCT{1},'r-');
hold on
plot(avrgACF_DATA_STRUCT{2},'g-');
hold on
plot(avrgACF_DATA_STRUCT{3},'m-');
hold on
plot(avrgACF_DATA_STRUCT{4},'b-');
hold on
grid on
yline(0)
xlabel("\tau_L (accepted MC steps)")
ylabel("Autocorrelation")
legend('\tau=0.08','\tau=0.8','\tau=0.4','\tau=40.0')
saveas(figure(1),"correlation-for-different-tau-linlin-plot.jpg")
%% LogLin Plot %%
figure(2)
semilogy(avrgACF_DATA_STRUCT{1},'ro-');
hold on
semilogy(avrgACF_DATA_STRUCT{2},'go-');
hold on
semilogy(avrgACF_DATA_STRUCT{3},'mo-');
hold on
semilogy(avrgACF_DATA_STRUCT{4},'bo-');
hold on
grid on
yline(0)
ylim([0.8, 1.0])
xlabel("\tau_L (accepted MC steps)")
ylabel("Autocorrelation")
hold on
%% Correlation time calculation and fitting %%
xStart = 1;
xCutoff = [8, 11 , 11, 11];
P_coefficients = {};
for i=1 : length(FILENAMES)
    P_coefficients{i} = polyfit( (xStart:xCutoff(i)), avrgACF_DATA_STRUCT{i}(xStart:xCutoff(i)), 1);
    corrTime = -1/(P_coefficients{i}(1));
    disp("Correlation time for tau = " + taus(i) + " is " + corrTime);
end
xx = 1:1:numLags;
fit = P_coefficients{4}(1)*xx+P_coefficients{4}(2); 
semilogy(fit,'k');
legend('\tau=0.08','\tau=0.8','\tau=0.4','\tau=40.0',' ')
saveas(figure(2),"correlation-for-different-tau-loglin-plot.jpg")
