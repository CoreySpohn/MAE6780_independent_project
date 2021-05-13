clear noises
%These are demonstrations
noise=make_noise(100,100)
get_noise(20,20,noise,100,100,1)

for i=1:500
    for j=1:500
        noises(i,j) = get_noise(i,j,noise,500,500,20);
    end
end
figure
imagesc(noises); colormap gray;

function val = get_noise(x,y,noise,max_x,max_y,scale_factor)

    %Input desired x and y values.
    %Input noise map.
    %Input largest x and y values.
    %95th percentile is such that the average of the 2.5 and 97.5 percentiles
    %are between 1 and -1, use scalefactor to change to appropriate scale.

    ptile=prctile(noise,[2.5,97.5],'all');
    scale_to_1=(abs(ptile(1))+abs(ptile(2)))/2;
    [n, m] = size(noise);
    y_s=n/max_y*y;
    x_s=m/max_x*x;
    val = interp2(noise,x_s,y_s);
    val = val*scale_factor/scale_to_1;

end


function noise = make_noise(n,m)
    %Input grid size, perlin noise map will be output.
    noise = zeros(n, m);
    noise = perlin_noise(noise);
end

function im = perlin_noise(im)
    %Generates smooth noise
    [n, m] = size(im);
    i = 0;
    w = sqrt(n*m);

    while w > 3
        i = i + 1;
        d = interp2(randn(n, m), i-1, 'spline');
        im = im + i * d(1:n, 1:m);
        w = w - ceil(w/2 - 1);
    end
end 