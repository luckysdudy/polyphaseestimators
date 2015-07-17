allocatemem(2^26)

N = 199
m = 3

p = m+1
V = matrix(N,p,i,j,i^(j-1))
K = mattranspose(V)*V
Ki = matsolve(K,matid(p))
Q = matid(N) - V*Ki*mattranspose(V)
print("Computing LLL")
T = qflllgram(Q)
R = Q*T

print("Lattice dimension is "matsize(R)[2])

basisfilename = Str("LSU_basis_N_"N, "_m_", m)
system("rm -f "basisfilename);
print("Writting basis to file "basisfilename)
write(basisfilename,R*1.0)

transformfilename = Str("LSU_transform_N_", N, "_m_", m)
system("rm -f "transformfilename);
print("Writting transformation to file "transformfilename)
write(transformfilename,T)

quit
