
####################
# Test eQTL Recipe #
####################

# load data 
bulklmmdir = dirname(pathof(BulkLMM));

pheno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-pheno-nomissing.csv");
pheno = BulkLMM.DelimitedFiles.readdlm(pheno_file, ',', header = false);
pheno_processed = pheno[2:end, 2:(end-1)].*1.0; # exclude the header, the first (transcript ID)and the last columns (sex)

geno_file = joinpath(bulklmmdir,"..", "data", "bxdData", "spleen-bxd-genoprob.csv")
geno = BulkLMM.DelimitedFiles.readdlm(geno_file, ',', header = false);
geno_processed = geno[2:end, 1:2:end] .* 1.0;

gmap_file = joinpath(bulklmmdir, "..", "data", "bxdData", "gmap.csv");
gInfo = BulkLMM.CSV.read(gmap_file, BulkLMM.DataFrames.DataFrame);

phenocovar_file = joinpath(bulklmmdir,"..", "data", "bxdData","phenocovar.csv");
pInfo = BulkLMM.CSV.read(phenocovar_file, BulkLMM.DataFrames.DataFrame);

results_path = joinpath(@__DIR__, "..", "data", "bxd", "multipletraits_results.he")
multipletraits_results = Helium.readhe(results_path);

kinship = calcKinship(geno_processed);

# use get_eQTL_accMb to get eQTL plotting inputs
x, y, z, mysteps, mychr = BigRiverPlots.get_eQTL_accMb(multipletraits_results, pInfo, gInfo; threshold = 5.0);

# generate plotting and save image as png to compare with the reference image 
plot_eQTL(multipletraits_results, pInfo, gInfo; threshold = 5.0)
savefig(joinpath(@__DIR__, "eQTL_new.png") )

img_test = FileIO.load(joinpath(@__DIR__,"..", "images","eQTL_test.png")); # ref image
img_new = FileIO.load(joinpath(@__DIR__, "eQTL_new.png")); # new image

# test plotting results
println("eQTL plot image test: ", @test img_test == img_new);

# clear new plot
rm(joinpath(@__DIR__, "eQTL_new.png"))

# testing plotting attributes
plot_obj = eqtlplot(x, y, z, mysteps, mychr);
println("eQTL plot attributes :x test: ", @test plot_obj[1][3].plotattributes[:x] == x);
println("eQTL plot attributes :y test: ", @test plot_obj[1][3].plotattributes[:y] == y);
println("eQTL plot attributes :z test: ", @test plot_obj[1][3].plotattributes[:marker_z] == z);


###################
# Test QTL Recipe #
###################

# Preprocessing 
traitID = 1112;
pheno_y = reshape(pheno_processed[:, traitID], :, 1);

# Scan 
single_results_perms = scan(
                        pheno_y, 
                        geno_processed, 
                        kinship; 
                        permutation_test = true, 
                        nperms = 2000, 
                        original = false
);                          


x, y, vecSteps, v_chr_names = get_plot_QTL_inputs(single_results_perms, gInfo)


# generate plotting and save image as png to compare with the reference image 
plot_QTL(single_results_perms, gInfo);
savefig(joinpath(@__DIR__, "QTL_new.png"));

plot_QTL(single_results_perms, gInfo, thresholds= [0.05, 0.75]);
savefig(joinpath(@__DIR__, "QTL_thrs_new.png"));


img_test = FileIO.load(joinpath(@__DIR__,"..", "images","QTL_test.png")); # ref image
img_thrs_test = FileIO.load(joinpath(@__DIR__,"..", "images","QTL_thrs_test.png")); # ref image with thresholds

img_new = FileIO.load(joinpath(@__DIR__, "QTL_new.png")); # new image
img_thrs_new = FileIO.load(joinpath(@__DIR__, "QTL_thrs_new.png")); # new image with thresholds

# test plotting results
println("QTL plot image test: ", @test img_test == img_new);
println("QTL plot image with thresholds test: ", @test img_thrs_test == img_thrs_new);

# clear new plot
rm(joinpath(@__DIR__, "QTL_new.png"))
rm(joinpath(@__DIR__, "QTL_thrs_new.png"))

# testing plotting attributes
max_lods = vec(mapslices(x -> maximum(x), single_results_perms; dims = 1));
thrs = map(x -> quantile(max_lods, x), [0.05, 0.75]);
plot_obj = qtlplot(x,y, vecSteps, v_chr_names, thrs);

idx_not_Inf = findall(x .!=Inf);
println("QTL plot attributes :x test: ", @test plot_obj[1][3].plotattributes[:x][idx_not_Inf] == x[idx_not_Inf]);
idx_not_Inf = findall(y .!=Inf);
println("QTL plot attributes :y test: ", @test plot_obj[1][3].plotattributes[:y][idx_not_Inf] == y[idx_not_Inf]);
