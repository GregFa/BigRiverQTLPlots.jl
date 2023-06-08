

####################
# Test eQTL Recipe #
####################

# load data 
bulklmmdir = dirname(pathof(BulkLMM));

gmap_file = joinpath(bulklmmdir, "..", "data", "bxdData", "gmap.csv");
gInfo = BulkLMM.CSV.read(gmap_file, BulkLMM.DataFrames.DataFrame);

# idx_geno = findall(occursin.(gInfo.Chr, "1 2 3 4 5"));
# gInfo_subset = gInfo[idx_geno, :];

phenocovar_file = joinpath(bulklmmdir, "..", "data", "bxdData", "phenocovar.csv");
pInfo = BulkLMM.CSV.read(phenocovar_file, BulkLMM.DataFrames.DataFrame);
# idx_pheno = findall(occursin.(pInfo.Chr, "1 2 3 4 5"));
# pInfo_subset = pInfo[idx_pheno, :];

pheno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-pheno-nomissing.csv");
pheno = BulkLMM.DelimitedFiles.readdlm(pheno_file, ',', header = false);
pheno_processed = pheno[2:end, 2:(end-1)] .* 1.0; # exclude the header, the first (transcript ID)and the last columns (sex)

# pheno_processed_subset = pheno_processed[:, idx_pheno];

geno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-bxd-genoprob.csv")
geno = BulkLMM.DelimitedFiles.readdlm(geno_file, ',', header = false);
geno_processed = geno[2:end, 1:2:end] .* 1.0;
# geno_processed_subset = geno_processed[:, idx_geno];

# kinship_subset = calcKinship(geno_processed_subset);


# multipletraits_results, heritability_results = bulkscan_null(
# 	pheno_processed_subset,
# 	geno_processed_subset,
# 	kinship_subset,
# );

# # use get_eQTL_accMb to get eQTL plotting inputs
# x, y, z, mysteps, mychr = BigRiverQTLPlots.get_eQTL_accMb(
# 	multipletraits_results,
# 	pInfo_subset,
# 	gInfo_subset;
# 	threshold = 5.0,
# );

# # generate plotting and save image as png to compare with the reference image 
# plot_eQTL(multipletraits_results, pInfo_subset, gInfo_subset; threshold = 5.0)
# savefig(joinpath(@__DIR__, "eQTL_new.png"))

# img_test = FileIO.load(joinpath(@__DIR__, "..", "images", "eQTL_test.png")); # ref image
# img_new = FileIO.load(joinpath(@__DIR__, "eQTL_new.png")); # new image

# # test plotting results
# println("eQTL plot image test: ", @test img_test == img_new);

# # clear new plot
# # rm(joinpath(@__DIR__, "eQTL_new.png"))

# # testing plotting attributes
# plot_obj = eqtlplot(x, y, z, mysteps, mychr);
# println("eQTL plot attributes :x test: ", @test plot_obj[1][3].plotattributes[:x] == x);
# println("eQTL plot attributes :y test: ", @test plot_obj[1][3].plotattributes[:y] == y);
# println("eQTL plot attributes :z test: ", @test plot_obj[1][3].plotattributes[:marker_z] == z);


###################
# Test QTL Recipe #
###################

# Preprocessing 
traitID = 1112;
pheno_y = reshape(pheno_processed[:, traitID], :, 1);

kinship = calcKinship(geno_processed);
Helium.writehe(kinship, joinpath(@__DIR__, "K_test.he"))


# Scan 
single_results_perms = scan(
	pheno_y,
	geno_processed,
	kinship;
	permutation_test = true,
	nperms = 2000,
	original = true,
);

Helium.writehe(single_results_perms, joinpath(@__DIR__, "scan_perms.he"))

single_results = scan(
	pheno_y,
	geno_processed,
	kinship,
);

Helium.writehe(reshape(single_results.lod, :,1), joinpath(@__DIR__, "scan_test.he"))


x = [1,2,3,4,5,6,7,8,9,10].*1.0;
mX = reshape(x, :,1)
mX_perms = BulkLMM.transform_permute(mX; nperms = 2000, original = false);
Helium.writehe(mX_perms, joinpath(@__DIR__, "mX_perms.he"))

K_eigen = BulkLMM.eigen(kinship);
K_eigen = reshape(K_eigen.values,:,1);
Helium.writehe(K_eigen, joinpath(@__DIR__, "eigen_test.he"))


x, y, vecSteps, v_chr_names = get_plot_QTL_inputs(single_results.lod, gInfo);


# generate plotting and save image as png to compare with the reference image 
plot_QTL(single_results.lod, gInfo);
savefig(joinpath(@__DIR__, "QTL_new.png"));


thrs = BigRiverQTLPlots.perms_thresholds(single_results_perms, [0.90, 0.95])
plot_QTL(single_results.lod, gInfo, thresholds = thrs);
savefig(joinpath(@__DIR__, "QTL_thrs_new.png"));


img_test = FileIO.load(joinpath(@__DIR__, "..", "images", "QTL_test.png")); # ref image
img_thrs_test = FileIO.load(joinpath(@__DIR__, "..", "images", "QTL_thrs_test.png")); # ref image with thresholds

img_new = FileIO.load(joinpath(@__DIR__, "QTL_new.png")); # new image
img_thrs_new = FileIO.load(joinpath(@__DIR__, "QTL_thrs_new.png")); # new image with thresholds

# test plotting results
println("QTL plot image test: ", @test (img_test == img_new));
# println("QTL plot image test2: ", @test (sum(img_test .== img_new) == size(img_new,1)*size(img_new,2)));
# println("QTL plot image with thresholds test: ", @test img_thrs_test == img_thrs_new);
# println("QTL plot image test2: ", 
# 	@test (sum(img_thrs_test .== img_thrs_new) == size(img_thrs_new,1)*size(img_thrs_new,2)));

# clear new plot
# rm(joinpath(@__DIR__, "QTL_new.png"))
# rm(joinpath(@__DIR__, "QTL_thrs_new.png"))

# testing plotting attributes
plot_obj = qtlplot(x, y, vecSteps, v_chr_names, thrs);

## test draw a number from test:
rng = MersenneTwister(0);
test_drawn = rand(rng);
println("Number drawn from test: ", test_drawn)

## test first element from test:
println("First element from test: ", single_results_perms[1, 1])



idx_not_Inf = findall(x .!= Inf);
println("QTL plot attributes :x test: ", @test plot_obj[1][3].plotattributes[:x][idx_not_Inf] == x[idx_not_Inf]);
idx_not_Inf = findall(y .!= Inf);
println("QTL plot attributes :y test: ", @test plot_obj[1][3].plotattributes[:y][idx_not_Inf] == y[idx_not_Inf]);
