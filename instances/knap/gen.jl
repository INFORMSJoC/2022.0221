using IterativePWLB
for n in [10, 20, 50, 100, 200, 500, 1000], k in 1 : 10
	IterativePWLB.write_knapsack_trindade(n, "/home/ccontard/git/ipwlb/code/instances/knap/CN23/cn23_$(n)_$(k).dat", k)
end
