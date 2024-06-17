load('/Users/matthew.baron/Perforce/audio-engineering/DSP/project_apollo/P1aFA4ArrayDesign/polars/ApolloP1a_FA4_20181031_raws.mat')
X1 = fft(x1);
figure; polarMap(sum(X1, 3), fs, ang(:,1), 1/12, [100 10000], [], 30, 2, 'none'); colormap('jet')
figure; polarMap(fft(sum(x1, 3)), fs, ang(:,1), 1/12, [100 10000], [], 30, 2, 'none'); colormap('jet')
figure; polarMap(fft(sum(y1, 3)), fs, ang(:,1), 1/12, [100 10000], [], 30, 2, 'none'); colormap('jet')
figure; polarMap(fft(sum(x1, 3)), fs, ang(:,1), 1/12, [100 10000], [], 30, 2, 'none'); colormap('jet')