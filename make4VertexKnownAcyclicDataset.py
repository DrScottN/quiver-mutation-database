from generation import *
M1=[[0, -4, 0, 1],
	[4, 0, 0, -2],
	[0, 0, 0, -2],
	[-1, 2, 2, 0]]

M2=[[0, -4, 0, 1],
	[4, 0, -1, 0],
	[0, 1, 0, -2],
	[-1, 0, 2, 0]]

M3=[[0, -3, 0, 1],
	[3, 0, 0, -1],
	[0, 0, 0, -2],
	[-1, 1, 2, 0]]

M4=[[0, -2, -1, 1],
	[2, 0, -1, 1],
	[1, 1, 0, -2],
	[-1, -1, 2, 0]]

M5=[[0, -2, -1, 1],
	[2, 0, 0, 0],
	[1, 0, 0, -2],
	[-1, 0, 2, 0]]

M6=[[0, -3, 0, 1],
	[3, 0, -1, 0],
	[0, 1, 0, -2],
	[-1, 0, 2, 0]]

M7=[[0, -2, -1, 2],
	[2, 0, -1, 0],
	[1, 1, 0, -2],
	[-2, 0, 2, 0]]

M8=[[0, -3, 0, 2],
	[3, 0, -2, -2],
	[0, 2, 0, -1],
	[-2, 2, 1, 0]]

M9 = [[0, -2, 0, 1],
	[2, 0, -1, 0],
	[0, 1, 0, -2],
	[-1, 0, 2, 0]]
provenCyclicIso = [boxQuiver(2,2), dreadedTorus(), Quiver(M1), Quiver(M2), Quiver(M3), Quiver(M4), Quiver(M5), Quiver(M6), Quiver(M7), Quiver(M8), Quiver(M9)]
for Q in provenCyclicIso:
    print(Q.matrix)
    print(Q.alexander_poly())
    print()

provenCyclic = []
for Q in provenCyclicIso:
    for p in permutations(4):
        Qp = isomorphicQuiver(Q, p)
        if Qp not in provenCyclic:
            provenCyclic.append(Qp)
provenCyclic = [mutationClass(Q) for Q in provenCyclic]
depth = 5
for QC in provenCyclic:
    for d in range(depth):
        QC.update()
    QC.mutationAcyclic=False #we proved this!
S = generateClassesWithKnownAcyclicity(4, 2, depth)
writeMutationClasses(itertools.chain(S, provenCyclic), "4vert2Max5MutationsAcyclicLabeled.csv", indexLabel=False)