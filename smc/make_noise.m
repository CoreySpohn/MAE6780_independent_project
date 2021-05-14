function noise = make_noise(n,m)
    %Input grid size, perlin noise map will be output.
    noise = zeros(n, m);
    noise = perlin_noise(noise);
end