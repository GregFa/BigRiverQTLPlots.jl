

#############
# Load Data #
#############

# load data 
bulklmmdir = dirname(pathof(BulkLMM));

gmap_file = joinpath(bulklmmdir, "..", "data", "bxdData", "gmap.csv");
gInfo = BulkLMM.CSV.read(gmap_file, BulkLMM.DataFrames.DataFrame);
idx_geno = findall(occursin.(gInfo.Chr, "1 2 3 4 5"));
gInfo_subset = gInfo[idx_geno, :];

phenocovar_file = joinpath(bulklmmdir, "..", "data", "bxdData", "phenocovar.csv");
pInfo = BulkLMM.CSV.read(phenocovar_file, BulkLMM.DataFrames.DataFrame);
idx_pheno = findall(occursin.(pInfo.Chr, "1 2 3 4 5"));
pInfo_subset = pInfo[idx_pheno, :];

pheno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-pheno-nomissing.csv");
pheno = BulkLMM.DelimitedFiles.readdlm(pheno_file, ',', header = false);
pheno_processed = pheno[2:end, 2:(end-1)] .* 1.0; # exclude the header, the first (transcript ID)and the last columns (sex)
pheno_processed_subset = pheno_processed[:, idx_pheno];

geno_file = joinpath(bulklmmdir, "..", "data", "bxdData", "spleen-bxd-genoprob.csv")
geno = BulkLMM.DelimitedFiles.readdlm(geno_file, ',', header = false);
geno_processed = geno[2:end, 1:2:end] .* 1.0;
geno_processed_subset = geno_processed[:, idx_geno];

####################
# Test eQTL Recipe #
####################

# Kinship 
kinship_subset = calcKinship(geno_processed_subset);

# Scan
multipletraits_results, heritability_results = bulkscan_null(
	pheno_processed_subset,
	geno_processed_subset,
	kinship_subset,
);

# use get_eQTL_accMb to get eQTL plotting inputs
x, y, z, mysteps, mychr = BigRiverQTLPlots.get_eQTL_accMb(
	multipletraits_results,
	pInfo_subset,
	gInfo_subset;
	threshold = 5.0,
);

# generate plotting and save image as png to compare with the reference image 
plot_eQTL(multipletraits_results, pInfo_subset, gInfo_subset; threshold = 5.0)
savefig(joinpath(@__DIR__, "eQTL_test.png"))

img_ref = FileIO.load(joinpath(@__DIR__, "..", "images", "eQTL_example.png")); # ref image
img_test = FileIO.load(joinpath(@__DIR__, "eQTL_test.png")); # new image

# test plotting results
println("eQTL plot image test: ", @test img_test == img_ref);

# clear new plot
rm(joinpath(@__DIR__, "eQTL_test.png"))

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

# Kinship 
kinship = calcKinship(geno_processed);
kinship = round.(kinship, digits = 12);

# Scan 
single_results_perms = scan(
	pheno_y,
	geno_processed,
	kinship;
	permutation_test = true,
	nperms = 1000,
);

single_results = scan(
	pheno_y,
	geno_processed,
	kinship,
);


x, y, vecSteps, v_chr_names = get_plot_QTL_inputs(single_results.lod, gInfo);

# generate plotting and save image as png to compare with the reference image 
plot_QTL(single_results.lod, gInfo);
savefig(joinpath(@__DIR__, "QTL_test.png"));

# generate plotting with thresholds obtain from perms_thresholds()
plot_QTL(single_results_perms, gInfo, thresholds = [0.90, 0.95]);
savefig(joinpath(@__DIR__, "QTL_thrs_test_1.png"));

# generate plotting with manual thresholds
thrs = BigRiverQTLPlots.perms_thresholds(single_results_perms.L_perms, [0.90, 0.95]);
println(thrs)
plot_QTL(single_results_perms.lod, gInfo, thresholds = thrs);
savefig(joinpath(@__DIR__, "QTL_thrs_test_2.png"));


img_ref = FileIO.load(joinpath(@__DIR__, "..", "images", "QTL_example.png")); # ref image
img_thrs_ref = FileIO.load(joinpath(@__DIR__, "..", "images", "QTL_thrs_example.png")); # ref image with thresholds

img_test = FileIO.load(joinpath(@__DIR__, "QTL_test.png")); # new image
img_thrs_test_1 = FileIO.load(joinpath(@__DIR__, "QTL_thrs_test_1.png")); # new image with thresholds
img_thrs_test_2 = FileIO.load(joinpath(@__DIR__, "QTL_thrs_test_2.png")); # new image with thresholds

# test plotting results
println("QTL plot image test: ", @test (img_test == img_ref));
println("QTL plot image with thresholds (auto) test: ", 
@test sum(1 .*(img_thrs_test_1 .== img_thrs_ref))==size(img_thrs_ref,1)*size(img_thrs_ref,2));
println("QTL plot image with thresholds (auto) test: ", @test img_thrs_test_1 == img_thrs_ref);
println("QTL plot image with thresholds (manual vs auto) test: ", @test img_thrs_test_2 == img_thrs_test_1);

# clear new plot
rm(joinpath(@__DIR__, "QTL_test.png"))
# rm(joinpath(@__DIR__, "QTL_thrs_test_1.png"))
# rm(joinpath(@__DIR__, "QTL_thrs_test_2.png"))


# testing plotting attributes
plot_obj = qtlplot(x, y, vecSteps, v_chr_names, thrs);

idx_not_Inf = findall(x .!= Inf);
println("QTL plot attributes :x test: ", 
			@test plot_obj[1][3].plotattributes[:x][idx_not_Inf] == x[idx_not_Inf]);
idx_not_Inf = findall(y .!= Inf);
println("QTL plot attributes :y test: ", 
			@test plot_obj[1][3].plotattributes[:y][idx_not_Inf] == y[idx_not_Inf]);
