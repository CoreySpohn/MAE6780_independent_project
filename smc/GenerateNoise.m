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