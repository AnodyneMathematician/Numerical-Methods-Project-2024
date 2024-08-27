#= 
Here is the initial/deprecated version of the Arnoldi Algorithm.
A second version will be uploaded shortly which includes tolerance.

=#

# Nota Bene: In the Second Version A Number of Code Smells still exist, and will need to be fixed in Code Refactorization 
# (This includes parts of code which includes duplicated functionality)
using LinearAlgebra


#=
Paramaters
A - matrix whose eigenvalues we are investigating
B - random normalized vector
k - size of the Arnoldi Block

Outputs a Hessenberg Matrix via Arnoldi.
This Hessenberg Matrix outputs Ritz Eigenvalues which are obtained using
QR methods (this took me way too long to find out).
=#

#See bottom of notebook for currently debugging version which also includes
#Tolerance

# n - Matrix size 

function normalizedrandomvector(n,k)
    b = zeros(n,k)
    b[:,1] = rand(n)
    b[:,1] /= norm(b[:,1])
    return b
end


function arnoldihessenberg(A,b,k)
    #Initializing Matrices
    H = zeros(k, k)
    r = copy(b[:,1])
    b[:,1] = r / norm(r)
    # Arnoldi Iteration Process
    for j = 1:k
        if j > 1
            b[:,j] = r / H[j,j-1]
        end
        r = A * b[:,j] 
        # Performing Modified Gram Schmidt
        for i=1:j
            H[i,j] = dot(b[:,i], r)
            r -= H[i,j] * b[:,i]
        end

        if j<k
            H[j+1,j] = norm(r)
        end
    end
    return H
end

