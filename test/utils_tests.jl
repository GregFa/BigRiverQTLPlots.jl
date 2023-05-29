
####################
# TEST pseudoticks #
####################

vticks = sort(rand(10).*10) 
vticks = [1,4,5,7] |> float
pseudoticks(vticks)

println("Get pseudoticks test: ", @test pseudoticks(vticks) == [0.5, 2.5, 4.5, 6]);

####################
# TEST sortnatural #
####################

myvec = ["12", "2", "1", "Y", "10", "X"];

println("Sort natural test: ", @test sortnatural(myvec) == ["1", "2", "10", "12", "X", "Y"]);


