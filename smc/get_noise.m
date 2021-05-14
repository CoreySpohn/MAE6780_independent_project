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