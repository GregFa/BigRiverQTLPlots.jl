using RCall, Statistics, DelimitedFiles, Distributed,LinearAlgebra # StatsBase, 
using Downloads
setSeed(123)

#####################
# Gough Island Mice #
#####################

########
# Data #
########


markerinfo2d_file = joinpath("example", "goughF2_markerinfo_2d_step0.csv");
markerinfo2d = readdlm(markerinfo_file,',';skipstart=1);

idx_2d=findall(markerinfo2d[:,2].!="X");
mar2d=markerinfo2d[idx_2d,:]



# dropping markers whose distance is less than 0.25 (filtering)
incl_idx=findall(diff(mar2d[:,3]).>=0.25)




gInfo = gInfo = DataFrame(
    Locus = markerinfo2d[:, 1];
    Chr = string.(markerinfo2d[:, 2]), 
    Mb = markerinfo2d[:,3], 
);



# increase timeout response from the server
downloader = Downloads.Downloader()
downloader.easy_hook = (easy, info) -> Downloads.Curl.setopt(easy, Downloads.Curl.CURLOPT_LOW_SPEED_TIME, 300)

lod_url = "https://github.com/senresearch/FlxQTL-SI/raw/main/Result/goughf2_2dlod_chr7.txt"
lod_file = Downloads.download(lod_url, tempname(); downloader = downloader)
ch7 = readdlm(lod_file);

lod_url = "https://github.com/senresearch/FlxQTL-SI/raw/main/Result/goughf2_2dlod_chr8.txt"
lod_file = Downloads.download(lod_url, tempname(); downloader = downloader)
ch8 = readdlm(lod_file);

lod_url = "https://github.com/senresearch/FlxQTL-SI/raw/main/Result/goughf2_2dlod_chr10.txt"
lod_file = Downloads.download(lod_url, tempname(); downloader = downloader)
ch10 = readdlm(lod_file);