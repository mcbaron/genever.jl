# matPolarPlot.jl
# Filter and call into matlab for polar plots

using MATLAB, MAT

function matPolarPlot(x1, eqFIR, fs, ang)

    y1 = Array{Float64}(undef, 8192, 72, 8)
    for i in 1:numSpkr
        y1[:,:,i] = filt(db2amp(23) * eqFIR[:,i], [1], x1[:,:,1,1,i])
    end

    mat"""
    $f = (0:(size($x1,1)/2-1))/(size($x1,1)/2)*$fs/2;
    X1 = fft($y1);
    drv = [1:8];
    drv_axis = repmat([0 360 250 5500],9,1);
    figure
    for i=1:length(drv)
        subplot(3,3,i);
        if 0
            polarMap(X1(:,:,drv(i)),$fs,$ang(:,1),1/12,[100 10000], [],30,2, 'none')
            colormap('jet');
        else
            h = surf($ang(:,1),$f,20*log10(squeeze(abs(X1(1:(size($x1,1)/2),:,drv(i))))));
            shading interp
            axis(drv_axis(i,:));
            view(90,90)
            set(gca,'YScale','log');
            colormap('jet');
            caxis([-30 0]-20);
        end
        title(sprintf('M%d',i));
    end
    """

    # write out to matfile for s&g
    file = matopen("ApolloLego_001.mat", "w")
    write(file, "y1", y1)
    close(file)
end
