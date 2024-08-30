#Computes Line Graphs, and Contour Plots
using PlotlyJS
using BenchmarkTools

#Localized to my machine, will need to change these paths for this code to work on your machines.
include("C:/University/Random Matrix 26 08 2024.jl")
include("C:/University/QRHessFrancis 26 08 2024.jl")
include("C:/University/Arnoldi 26 08 2024.jl")


function hessfrancisqralg2(M::Matrix{Float64})
    
    #count = 0
    M = ComplexF64.(M) #Converts matrix into complex form
    n = size(M, 1) #Calcaulates size of matrix (assumed square)
    ϵ1 = 1e-3 #Error term (used to check when algorithm should terminate)
    ϵ2 = (ϵ1)^2 #Sqaure of Error term (used to avoid having to take square root of complex numbers)
    #Reduces matrix to hessenberg form
    newA = francisqrdecomp(M, n, ϵ1) #Performs francis double shift QR decomposition
    Eigen = extracteigen(M, n, ϵ1) #Makes list of eigenvalues of A
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


#Plots Matrix Size v. Sparsity

#Computes Meantimes/Standard Deviation for various levels of sparsity.
meantime = []
meantimesample = []
sparsemeantime = []
standarddeviation = []
#No of Samples for a given Sparsity and Matrix Size
ksample = 1
#Maximum Matrix Size
msize = 70
#Maximum Sparsity
nsparse = 40
#Step change in Sparsity
s1 = 5
#Step change in Matrix Size
s2 = 5

meantime = []
for sizes = 1:(Int64(msize/s2))
    #Determines size of a matrix for a step/iteration in the for loop
    stepsize = s2*sizes
    #Re-initializes sparse list for each change in matrix size
    sparsemeantime = []
    for sparsity = 1:Int(nsparse/s1)
        #Re-initializes sample list for each change in sparsity
        samplemeantime = []
        for samples = 1:ksample
            #Records time taken to compute a sample
            t1 = time()
            #Computes Arnoldi Hessenberg Function with Hess Francis/QR at various levels of percentage sparsity from s1% to nsparse%
            
            #Computes Arnoldi Hessenberg Function with Hess Francis/QR for different samples at various levels of percentage sparsity from s1% to nsparse%
            #H = hessfrancisqralg2(arnoldihessenberg(randmatrix(stepsize,Int64(floor(s1*sparsity*(stepsize)/100))),normalizedrandomvector(stepsize,Int64(ceil(stepsize/2))),Int64(ceil(stepsize/2))))
            #Computes Hess Francis/QR for different samples at various levels of percentage sparsity from s1% to nsparse%
            H = hessfrancisqralg2(randmatrix(stepsize,Int64(floor(s1*sparsity*(stepsize)/100))))
            elapsed_time = time() - t1;
            #Adds recorded value to list of samples
            #Removes 0 from median scores
            if elapsed_time > 0
                push!(samplemeantime,Float64(elapsed_time/stepsize^2))
            end

        end
        #Provides the value of a pixel on the contour plot - by calculating from the samples (determined by ksample) the median time taken to compute. 
        #Adds to sparsity vector.
        #Adds a 0 if no other values record
        if samplemeantime  == []
            push!(samplemeantime,0.0)
        end
        push!(sparsemeantime,median(samplemeantime))
        println(meantime)
    end
    #Provides the values of a line of a contour plot - All sparsities for a given Matrix Size.
    #Adds to meantime matrix.
    push!(meantime,Float64.(sparsemeantime))
end

A = meantime
#Combines All of the vectors of meantime into one singular matrix
a = []
a = Float64.(vcat(a, A[1]))
#Each Index Note
for i in 2:length(A)
    a = hcat(a,A[i])
end 
#For Logarithmic Plots
#=
loga = log.(a)
b = loga[1:10,1:20]
=#
a
#Plots subset of Matrix as Domain
#Change these values for what subset of the data you'd like to look at.
b = a[1:8,1:10]
#Produces Contour Plot
function contour4()
    z = b
    data = contour(;z=z, colorscale="Jet", connectgaps = true)

    layout = Layout(;title="Contour Plot for Sparsity v. Matrix Size - QR",xaxis=attr(title="Matrix Size"),yaxis=attr(title="Sparsity"))
    plot(data, layout)
end
contour4()

a2
a3 = a























a2 = a

a1

#Noisy Regime 0-40
#Samples 3, Logarithmic, Sparsity 0-100%, Matrix Size 5-600

#Plots No. of Eigenvalues v. Sparsity