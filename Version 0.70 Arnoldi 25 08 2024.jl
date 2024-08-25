#= 
Here is the initial/deprecated version of the Arnoldi Algorithm.
A second version will be uploaded shortly which includes tolerance.

=#

# Nota Bene: In the Second Version A Number of Code Smells still exist, and will need to be fixed in Code Refactorization 
# (This includes parts of code which includes duplicated functionality)


using Plots
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



#Temporary Function to calculate Eigenvalues before linked with Robbie's function

function eig_k(H,k)
    Hk = H[1:k,1:k]
    L = eigen(Hk).values # Eigenvalue estimates at step k        
    return L    
end



#Test Case 1
#Dimension of (Square) Matrix
n = 100
#Size of Arnoldi Block
k = n
#Initialization
X = rand(n,n)
# Random (or guessed) initial vector
b = zeros(n,k)
b[:,1] = rand(n)
b[:,1] /= norm(b[:,1])

# Arnoldi to Hessenberg
H = arnoldihessenberg(X,b,k)'
# Hessenberg to Eigenvalues (Call Robbie's function Here, in lieu:
# Temporary function is used)
Eig = eig_k(X,100)
#Plotting Eigenvalues as a scatter plot:
scatter(real(Eig),imag(Eig))
 



#Test Case 2

#To be copied from previous workbook.


#More Test Cases Should be Completed.


# A more readable verson, with greater functionality but with significant debugging issues:

#=
function arnoldi(A,b,k,tol)
    #Nomralization of b
    r = copy(b)
    #Initializing Matrices Q & H
    Q = zeros(length(r),k+1)
    println(Q)
    H = zeros(k+1,k)
    println(H)
    Q[:,1] = r
    Q[:,1] /= norm(Q[:,1])
    #Performing Iteration
    for j = 1:k
        Q[:,j+1] = A*Q[:,j]
        #Modified Gram Schmidt (use Robbies function?)
        for i=1:j
            H[i,k] = dot(Q[:,i], Q[:,j+1])
            Q[:,j+1] -= H[i,j] * Q[:,i]
            #Q[:,j+1] -= H[i,j] * Q[:,i]
        end
        H[j+1,j] = norm(Q[:,j+1])


        if abs(H[j+1,j]) < tol
            return (H[:j,:j],Q[:,:j])
        end
        #Normalization
        Q[:,j+1] /= H[j+1,j]
    end
    return H
end


=#


#Old Native Implementation of QR methods (severly deprecated)


#=


#using Random
#rng = MersenneTwister(18);

#Performs Givens Rotation on two columns of a matrix, to derive an Upper Triangular Matrix
function givensrotation(v1,v2)
    if v2 == 0
        c = 1
        s = 0
    else
        if abs(v2) > abs(v1)
            mu = -v1/v2
            s = 1.0/sqrt(1.0 + mu*mu)
            c = s*mu
        else
            mu = -v2/v1
            c = 1.0/sqrt(1.0 + mu*mu)
            s = c*mu
        end
    end
    return (c,s)
end

function hessenbergtriangular!(H,n)
    for k=1:n-1
        c, s = givensrotation(H[k,k], H[k+1,k])
        # Apply the Givens rotation to each row.
        for j=k:n
            H[k,j], H[k+1,j] =
                ( c * H[k,j] - s * H[k+1,j],
                s * H[k,j] + c * H[k+1,j] )
        end
    end
    return H
end

# Hessenberg to Triangular
C = hessenbergtriangular!(H,n)
#Triangular to Ritz
diag(C)

=#