#Please Note that file may be tricky to run in Jupyter Notebooks, 
#and may need to be run in .jl form in an IDE of your choice.

#Computes Line Graphs, and Contour Plots
using PlotlyJS
using BenchmarkTools
using ArnoldiMethod
using Statistics

#Localized to my machine, will need to change these paths for this code to work on your machines.
include("C:/University/Random Matrix 26 08 2024.jl")
include("C:/University/QRHessFrancis 26 08 2024.jl")
include("C:/University/Arnoldi 26 08 2024.jl")

#Modifed Function from the QRHessFranis to speed up computation of version 2
#If necessary.

function hessfrancisqralg2(M::Matrix{Float64})
    
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
        if count > 100
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
mediantime = []
meantimesample = []
sparsemediantime = []

sparsedeviation = []
deviation = []
#No of Samples for a given Sparsity and Matrix Size
ksample = 40
#Maximum Matrix Size
msize = 300
#Maximum Sparsity
nsparse = 90
#Step change in Sparsity
s1 = 5
#Step change in Matrix Size
s2 = 5

#Only relevant for Eigenvalue Versions of Plot 
#(Instead of Sparsity, or Matrix Size)
#Maximum No of Eigenvalues Computed (Portion of Total Matrix)
#eigenfraction = 1/4
#Change in Fraction of Eigenvalues Computed
#eigenchange = 1/10


mediantime = []
deviation = []
println("Plot Starting!")
for sizes = 1:(Int64(msize/s2))
    #Determines size of a matrix for a step/iteration in the for loop
    stepsize = s2*sizes
    #Re-initializes sparse list for each change in matrix size
    sparsedeviation = []
    sparsemediantime = []
    for step = 1:20
        #Re-initializes sample list for each change in sparsity
        samplemediantime = []
        for samples = 1:ksample
            #Records time taken to compute a sample 

            #Version 1: Speed of Arnoldi
            #t1 = time()
            #Computes Arnoldi Hessenberg Function for different samples of a given percentage sparsity from s1% to nsparse%
            #K = arnoldihessenberg(randmatrix(stepsize,Int64(floor(s1*30*(stepsize)/100))),normalizedrandomvector(stepsize,Int64(ceil(stepsize*(eigenfraction+eigenstep*eigenchange)))),Int64(ceil(stepsize*eigenfraction+eigenstep*eigenchange)))
            #elapsed_time = time() - t1;

            #Version 2: Difference Between Time of Arnoldi vs. QR to
            #to compute a result. 
            
            #Computes Arnoldi Hessenberg Function with Hess Francis/QR for different samples at various levels of percentage sparsity from s1% to nsparse%
            #t1 = time()
            #H = hessfrancisqralg2(arnoldihessenberg(randmatrix(stepsize,Int64(floor(s1*step*(stepsize)/100))),normalizedrandomvector(stepsize,Int64(ceil(stepsize/2))),Int64(ceil(stepsize/2))))
            #elapsed_time = time() - t1;

            #Computes Hess Francis/QR for different samples at various levels of percentage sparsity from s1% to nsparse%
            #t2 = time()
            #H = hessfrancisqralg2(randmatrix(stepsize,Int64(floor(s1*step*(stepsize)/100))))
            #elapsed_time2 = time() - t2;

            #Version 3:
            #To show the functionality of the Contour Plots 
            #without making you spend an afternoon computing this
            #we have provided a simplifed version 2
            #using in house functions.
            #This version invloves comparing the Hessenberg, and Arnoldi 
            #Decompositions directly.
            #Should hopefully be done after a Cup of Coffee.
            A = randmatrix(stepsize,Int64(floor(s1*step*(stepsize)/100)))
            t1 = time()
            #Arnoldi Method Function - nev stands for number of eigenvalues
            #encapsulated with size(A,1) for when you have matrices 
            #with less than dimension 7.
            partialschur(A, nev = min(7,size(A,1)), tol=1e-2, which=:SR)
            elapsed_time = time() - t1;
            t2 = time()
            #Hessenberg Decomposition
            hessenberg(A)
            elapsed_time2 = time() - t2
            #Adds recorded value to list of samples
            elapsed_time = elapsed_time2 - elapsed_time
            #Removes 0 from median scores
            if elapsed_time != 0
                push!(samplemediantime,Float64(elapsed_time))
            end

        end
        #Provides the value of a pixel on the contour plot - by calculating from the samples (determined by ksample) the median time taken to compute. 
        #Adds to sparsity vector.
        #Adds a 0 if no other values recorded
        if samplemediantime  == []
            push!(samplemediantime,0.0)
        end
        #If you want to print something to make sure Julia's working/check progress.
        #println(size(sparsemediantime))
        #Computes Median of Samples Adds to Sparsity List for a given Matrix Size
        push!(sparsemediantime,median(samplemediantime))
        push!(sparsedeviation,std(samplemediantime))

    end
    #Provides the values of a line of a contour plot - All sparsities for a given Matrix Size.
    #Adds to mediantime matrix.

    push!(mediantime,Float64.(sparsemediantime))
    push!(deviation,Float64.(sparsedeviation))
end

println("Plot Finished!")
A = mediantime
B = deviation
#If you would like to map the standard deviation uncomment the below line

#A = B

#Combines All of the vectors of mediantime/deviation into one singular matrix
a = []
a = Float64.(vcat(a, A[1]))
#Each Index Note
for i in 2:length(A)
    a = hcat(a,A[i])
end 
#For Logarithmic Plots
#Uncomment the following (Please Note won't work for difference calculations 
#as we have negative values)

#loga = log.(a)
#b = loga[1:20,1:40]

#Checking No Untoward Behvaiour With Plot by Printing Out Matrix

a

#Plots subset of Matrix as Domain
#Change these values for what subset of the data you'd like to look at.
b = a[1:20,1:60]
#Use this incase there are any anomalous results which 
#you temporarily want to remove to get better insight into 
#the rest of the graph
#(usually from low samples/odd behaviour at index 1,1)
#b[1,1] = 0.000001
#Produces Contour Plot
function contour4()
    z = b
    data = contour(;z=z, colorscale="Jet", connectgaps = true, dy=1/20, y0=0.05, dx=5, x0=5,colorbar=attr(;title="Time (seconds)",titleside="right"))

    layout = Layout(;title="Contour Plot for Sparsity v. Matrix Size",xaxis=attr(title="Matrix Size"),yaxis=attr(title="Sparsity Proportion"))
    plot(data, layout)
end
contour4()
























a2 = a

a1

#Noisy Regime 0-40
#Samples 3, Logarithmic, Sparsity 0-100%, Matrix Size 5-600

#Plots No. of Eigenvalues v. Sparsity