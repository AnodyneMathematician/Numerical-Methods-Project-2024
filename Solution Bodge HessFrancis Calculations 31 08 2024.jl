#Locally imported. Paths will very likely be different for your machine.
include("C:/University/Random Matrix 26 08 2024.jl")
include("C:/University/QRHESSFRANCIS2 26 08 2024.jl")
include("C:/University/Arnoldi 26 08 2024.jl")

#Additional Count to 
errored = 0

#Please note the count in the code

#Varied Form of the Hess Francis Alg for termination when iteration is at large values.
function hessfrancisqralg2(M::Matrix{Float64})
    
    #count = 0
    M = ComplexF64.(M) #Converts matrix into complex form
    n = size(M, 1) #Calcaulates size of matrix (assumed square)
    ϵ1 = 1e-3 #Error term (used to check when algorithm should terminate)
    ϵ2 = (ϵ1)^2 #Sqaure of Error term (used to avoid having to take square root of complex numbers)
    #Reduces matrix to hessenberg form
    A = hessen(M, n, ϵ1) #Reduces matrix to hessenberg form

    newA = francisqrdecomp(A, n, ϵ1) #Performs francis double shift QR decomposition
    Eigen = extracteigen(A, n, ϵ1) #Makes list of eigenvalues of A
    newEigen = extracteigen(newA, n, ϵ1) #Makes list of eigenvalues of newA
    count = 0
    while all(abs2.(Eigen-newEigen) .< ϵ2) == false #checks if the eigenvalues of the newA and A are all within ϵ1
        count += 1
        #of each other and continues iterating if they are not
        A = newA #Updates A
        newA = francisqrdecomp(A, n, ϵ1) #Computes newA
        Eigen = newEigen #Updates list of eigenvalues of A
        newEigen = extracteigen(newA, n, ϵ1) #Updates list of eigenvalues of newA
        if count > 200
            global errored += 1
            break
        end
    end
    println(count)
    for i in 1:n #Iterates over all eigenvalues
        if abs2(newEigen[i].im) < ϵ2 #Checks if imaginary part of eigenvalue is within ϵ1 of 0
            newEigen[i] = ComplexF64(newEigen[i].re) #If so, removes imaginary component so that the eigenvalue is real
        end
        if abs2(newEigen[i].re) < ϵ2 #Checks if real part of eigenvale is within ϵ1 of 0
            newEigen[i] = ComplexF64(newEigen[i].im) #If so, removes real component so that the eigenvalue is imaginary
        end
    end
    return newEigen #Returns eigenvalues
end

count2 = 0

#Runs the algorithm above 100 times.
for i = 1:100
    count2 += 1
    #Uncomment this print statement if you want to see what iteration the for loop is computing.
    #println(count2)
    #Inside function from Arnoldi should be imported in-order to work.
    hessfrancisqralg2(arnoldihessenberg(randmatrix(40,40),normalizedrandomvector(40,40),40))
end
#Prints No of Matrices which failed to terminate within 200 iterations.
errored



