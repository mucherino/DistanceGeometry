
# The Stick Model for Distance Geometry in 2D
#
# when stick orientations are fixed, 
# the problem is convex quadratic with linear constraints
#
# April 19, 2025
#
# AM

using LinearAlgebra
using JuMP
using Ipopt
using Plots
using Printf

# generate an n-vector of random pairs of floats between 0 and 1
function gen_coordinates(n::Int) :: Vector{Tuple{Float64,Float64}}
   if (n <= 0) throw(ArgumentError("The size n must be a positive integer")) end
   points = Vector{Tuple{Float64,Float64}}()
   for i = 1:n
      push!(points,(rand(),rand()))
   end
   return points
end

# verifying whether a given Float64 matrix contains only positive entries
function is_positive(matrix::Matrix{Float64})
   for k in eachindex(matrix)
      if (matrix[k] < 0.0) return false end
   end
   return true
end

# verifying whether a given matrix is squared
function is_squared(matrix::Matrix{T}) where T
   n = size(matrix,1)
   m = size(matrix,2)
   return m == n
end

# verifying whether a given matrix is symmetric
function is_symmetric(matrix::Matrix{T}) where T
   if (!is_squared(matrix)) return false end
   n = size(matrix,1)
   for i = 1:n
      for j = 1:i - 1
         if (matrix[i,j] != matrix[j,i]) return false end
      end
   end
   return true
end

# verifying whether a given matrix is antisymmetric
function is_antisymmetric(matrix::Matrix{T}) where T <: Number
   if (!is_squared(matrix)) return false end
   n = size(matrix,1)
   for i = 1:n
      for j = 1:i - 1
         if (matrix[i,j] != -matrix[j,i]) return false end
      end
   end
   return true
end

# verifying whether a given Float64 matrix has off-diagonal zeros
function has_offdiagonal_zeros(matrix::Matrix{Float64})
   n = size(matrix,1)
   m = size(matrix,2)
   for i = 1:n
      for j = 1:i - 1
         if (matrix[i,j] == 0.0) return true end
      end
      for j = i + 1:m
         if (matrix[i,j] == 0.0) return true end
      end
   end
   return false
end

# computing the ratio between 
# - the number of off-diagonal non-zero elements and 
# - the total elements of the input Float64 symmetric matrix
function density(matrix::Matrix{Float64})
   if (!is_symmetric(matrix)) throw(ArgumentError("The input matrix is not symmetric")) end
   n = size(matrix,1)
   count = 0
   for i = 1:n
      for j = 1:i - 1
         if (matrix[i,j] != 0.0) count = count + 2 end
      end
   end
   return count/(n^2)
end

# computing the distance matrix from a given point configuration
function distance_matrix(points::Vector{Tuple{Float64,Float64}})
   n = length(points)
   distances = zeros(n,n)
   for i = 1:n
      for j = 1:n
         distances[i,j] = sqrt( (points[j][1] - points[i][1])^2 + (points[j][2] - points[i][2])^2 )
      end
   end

   # the distances are ready
   return distances
end

# computing the orientation matrix from a given point configuration (with its pre-computed distance matrix)
function orientation_matrix(points::Vector{Tuple{Float64,Float64}},distances::Matrix{Float64})
   if (!is_positive(distances)) throw(ArgumentError("The input distance matrix contains negative entries")) end
   if (!is_symmetric(distances)) throw(ArgumentError("The input matrix is not symmetric (is supposed to be a distance matrix)")) end
   n = length(points)
   if (n != size(distances,1)) throw(ArgumentError("Number of input points does not match with distance matrix size")) end
   orientations = zeros(ComplexF64,n,n)

   # computing the orientations
   for i = 1:n
      for j = 1:n
         if distances[i,j] == 0.0
            orientations[i,j] = complex(0.0,0.0)  # special case
         else
            diff_x = points[j][1] - points[i][1]
            diff_y = points[j][2] - points[i][2]
            orientations[i,j] = complex(diff_x / distances[i,j],diff_y / distances[i,j])
         end
      end
   end

   # the orientations are ready
   return orientations
end

# reconstructing the original point configuration from the computed distance and orientation matrices
function reconstruct(distances::Matrix{Float64},orientations::Matrix{ComplexF64})
   if (!is_positive(distances)) throw(ArgumentError("The input distance matrix contains negative entries")) end
   if (!is_symmetric(distances)) throw(ArgumentError("The first input matrix is not symmetric (is supposed to be a distance matrix)")) end
   n = size(distances,1)
   if (!is_antisymmetric(orientations)) throw(ArgumentError("The second input matrix is not antisymmetric (is supposed to be an orientation matrix)")) end
   if (n != size(orientations,1)) throw(ArgumentError("The two input matrices differ in size")) end 
   if (has_offdiagonal_zeros(distances)) throw(ArgumentError("Not all distances are available in the distance matrix")) end
   points = Vector{Tuple{Float64,Float64}}()

   # the first point is arbitrarly positioned in the origin (0,0)
   push!(points,(0.0,0.0))

   # computing the coordinates of all other points
   for i = 2:n
      # since we have all distances and orientations, we can choose any reference point
      # -> here we always use the first point (positioned in (0,0)), but the code is general
      x = points[1][1] + distances[1,i]*real(orientations[1,i])
      y = points[1][2] + distances[1,i]*imag(orientations[1,i])
      push!(points,(x,y))
   end

   # the configuration is ready
   return points
end

# perturbing the values of a symmetric matrix having positive Float64 entries
function perturb!(matrix::Matrix{Float64},epsilon::Float64)
   if (!is_positive(matrix)) throw(ArgumentError("The input matrix contains negative entries")) end
   if (!is_symmetric(matrix)) throw(ArgumentError("The input matrix is not symmetric")) end
   if (epsilon <= 0.0) throw(ArgumentError("The epsilon variable must be strictly positive")) end
   n = size(matrix,1)
   for i = 1:n
      for j = 1:i - 1
         delta = rand()*epsilon - 0.5*epsilon
         matrix[i,j] = matrix[i,j] + delta
         if (matrix[i,j] < 0.0)  matrix[i,j] = 0.0 end
         matrix[j,i] = matrix[i,j]
      end
   end
end

# perturbing the angles related to the complex numbers forming an anti-symmetric matrix
function perturb!(matrix::Matrix{ComplexF64},epsilon::Float64)
   if (!is_antisymmetric(matrix)) throw(ArgumentError("The input matrix is not antisymmetric")) end
   if (epsilon <= 0.0) throw(ArgumentError("The epsilon error to introduce must be strictly positive")) end
   if (epsilon > 2.0*pi) throw(ArgumentError("The epsilon error cannot be larger than 2pi")) end
   n = size(matrix,1)
   for i = 1:n
      for j = 1:i - 1
         if i != j
            current_angle = angle(matrix[i,j])
            delta = rand()*epsilon - 0.5*epsilon
            new_angle = current_angle + delta
	    matrix[i,j] = abs(matrix[i,j])*cis(new_angle)  # cis(x) is equivalent to e^{ix}
            matrix[j,i] = -matrix[i,j]  # antisymmetric (different from Harmitian)
         end
      end
   end
end

# sparsifying a symmetric matrix of positive Float64 entries, together with the corresponding orientation matrix
# -> in the Float64 matrix, all entries larger than threshold are replaced by zeros
# -> the operation is NOT performed if the current entry is the only remaining one which is nonzero in corresponding row and column
# -> if (i,j) is replaced by zero in the distance matrix, then we do the same in the ComplexF64 matrix
function sparsify!(distances::Matrix{Float64},orientations::Matrix{ComplexF64},threshold::Float64)
   if (!is_positive(distances)) throw(ArgumentError("The input distance matrix contains negative entries")) end
   if (!is_symmetric(distances)) throw(ArgumentError("The input distance matrix is not symmetric")) end
   if (!is_antisymmetric(orientations)) throw(ArgumentError("The input orientation matrix is not antisymmetric")) end
   n = size(distances,1)
   if (n != size(orientations,1)) throw(ArgumentError("The two input matrices differ in size")) end
   if (threshold <= 0.0) throw(ArgumentError("The threshold must be strictly positive")) end
   for i = 1:n
      for j = 1:i - 1
         if distances[i,j] > threshold
            value = distances[i,j]
            distances[i,j] = 0.0
            distances[j,i] = 0.0
            if all(distances[i,k] == 0.0 for k in 1:n if k != j) && all(distances[k,j] == 0.0 for k in 1:n if k != i)
               distances[i,j] = value
               distances[j,i] = value
            else
               orientations[i,j] = complex(0.0,0.0)
               orientations[j,i] = orientations[i,j]
            end
         end
      end
   end
end

# generating and solving the Stick Model with JuMP (and Ipopt as predefined solver)
function stick_model(distances::Matrix{Float64},orientations::Matrix{ComplexF64})
   if (!is_positive(distances)) throw(ArgumentError("The input distance matrix contains negative entries")) end
   if (!is_symmetric(distances)) throw(ArgumentError("The input distance matrix is not symmetric")) end
   if (!is_antisymmetric(orientations)) throw(ArgumentError("The input orientation matrix is not antisymmetric")) end
   n = size(distances,1)
   if (n != size(orientations,1)) throw(ArgumentError("The two input matrices differ in size")) end

   # creating the JuMP model with Ipopt
   model = Model(Ipopt.Optimizer)
   set_optimizer_attribute(model,"print_level",1)  # only errors and warnings

   # defining variables (one "stick" per distance)
   @variable(model,cx[1:n,1:n])
   @variable(model,cy[1:n,1:n])

   # defining the constraints (sticks cannot change length)
   for i = 1:n
      for j = 1:n
         if i != j
            if distances[i,j] != 0.0 && orientations[i,j] != complex(0.0,0.0)
               theta = angle(orientations[i,j])  # getting the angle (theta) from the complex orientation
               @constraint(model,cx[i,j] == cx[j,i] + distances[i,j]*cos(theta))  # x-component
               @constraint(model,cy[i,j] == cy[j,i] + distances[i,j]*sin(theta))  # y-component
            end
         end
      end
   end

   # defining the objective function
   # -> minimizing the difference between "copies" of the same point
   obj = 0.0
   for i = 1:n
      for j = 1:n
         for k = 1:n
            if j != k
               obj = obj + (cx[i,k] - cx[i,j])^2 + (cy[i,k] - cy[i,j])^2
            end
         end
      end
   end
   @objective(model,Min,obj)

   # calling the solver
   optimize!(model)

   # creating a new point configuration from a given solution to the Stick Model
   points = Vector{Tuple{Float64,Float64}}()
   for i = 1:n
      x = 0.0
      y = 0.0
      for j = 1:n
         if i != j
            x = x + value(cx[i,j])
            y = y + value(cy[i,j])
         end
      end
      x = x / (n - 1)
      y = y / (n - 1)
      push!(points,(x,y))
   end

   # extracting the solution (coordinates) and compute the average maximum error
   avg_errors = zeros(n)
   max_error = nothing
   for i = 1:n
      for j = 1:n
         if i != j
            avg_errors[i] = avg_errors[i] + (points[i][1] - value(cx[i,j]))^2 + (points[i][2] - value(cy[i,j]))^2
         end
      end
      avg_errors[i] = sqrt(avg_errors[i] / (n - 1))
      if max_error == nothing || avg_errors[i] > max_error
         max_error = avg_errors[i]
      end
   end

   # returning the maximum error, and the computed point configuration
   return max_error, points
end

# computing a random orientation matrix
function random_orientations(n::Int64)
   if (n <= 0) throw(ArgumentError("The size of the orientation matrix must be positive")) end
   orientations = zeros(ComplexF64,n,n)
   for i = 1:n
      for j = 1:i - 1
         angle = 2*pi*rand()
         orientations[i,j] =  complex(cos(angle),sin(angle))
         orientations[j,i] = -orientations[i,j]
      end
   end

   # the random orientations are ready
   return orientations
end

# Stick Heuristic
function stick_heuristic(maxit::Int64,distances::Matrix{Float64},epsilon::Float64)
   if (maxit <= 1) throw(ArgumentError("The maximum number of iterations must be strictly larger than 1 ")) end
   if (!is_positive(distances)) throw(ArgumentError("The input distance matrix contains negative entries")) end
   if (!is_symmetric(distances)) throw(ArgumentError("The input distance matrix is not symmetric")) end
   if (epsilon <= 0.0) throw(ArgumentError("The tolerance epsilon must be strictly positive")) end

   # main initializations
   n = size(distances,1)
   current = nothing
   error = nothing
   orientations = random_orientations(n)
   converged = false
   iterations = 1

   # main loop
   while !converged && iterations < maxit

      # solving the current Stick Problem
      unused, current = stick_model(distances,orientations)

      # computing the distance matrix from current solution
      current_distances = distance_matrix(current)

      # computing the orientation matrix from the current solution
      current_orientations = orientation_matrix(current,current_distances)

      # convergence?
      error = orientation_error(current_orientations,orientations)
      println("#",iterations," stick model ",error)
      if (error < epsilon)
         converged = true
      else
         orientations = current_orientations
         iterations = iterations + 1
      end
   end

   # returning the current solution (with its error estimation)
   return error, current
end

# computing the average distance error (zero entries are ignored)
function distance_error(dist1::Matrix{Float64},dist2::Matrix{Float64})
   if (!is_positive(dist1)) throw(ArgumentError("The first input matrix contains negative entries")) end
   if (!is_positive(dist2)) throw(ArgumentError("The second input matrix contains negative entries")) end
   if (!is_symmetric(dist1)) throw(ArgumentError("The first input matrix is not symmetric")) end
   if (!is_symmetric(dist2)) throw(ArgumentError("The second input matrix is not symmetric")) end
   n = size(dist1,1)
   if (n != size(dist2,1)) throw(ArgumentError("The two input matrices differ in size")) end

   # comparing the two distance matrices
   avg = 0.0
   count = 0
   for i = 1:n
      for j = 1:i - 1
         if dist1[i,j] != 0.0 && dist2[i,j] != 0.0
            avg = avg + (dist2[i,j] - dist1[i,j])^2
            count = count + 1
         end
      end
   end
   if (count != 0)  avg = sqrt(avg / count) end

   # average error on the distances
   return avg
end

# computing the average orientation error (zero entries are ignored)
function orientation_error(orient1::Matrix{ComplexF64},orient2::Matrix{ComplexF64})
   if (!is_antisymmetric(orient1)) throw(ArgumentError("The first input orientation matrix is not antisymmetric")) end
   if (!is_antisymmetric(orient2)) throw(ArgumentError("The second input orientation matrix is not antisymmetric")) end
   n = size(orient1,1)
   if (n != size(orient2,1)) throw(ArgumentError("The two input matrices differ in size")) end

   # comparing the two matrices of complex numbers
   avg = 0.0
   count = 0
   czero = complex(0.0,0.0)
   for i = 1:n
      for j = 1:i - 1
         if orient1[i,j] != czero && orient2[i,j] != czero
            value = (angle(orient2[i,j]) - angle(orient1[i,j]))^2
            theta = 0.0
            if (value != 0.0)
               theta  = sqrt(value)
               thetap = abs(pi - theta)
               if (thetap < theta) theta = thetap end
            end
            avg = avg + theta
            count = count + 1
         end
      end
   end
   if (count != 0)  avg = sqrt(avg / count) end

   # avergage error on the orientations
   return avg
end

# centering an input vector of point coordinates
# -> the modified points are returned in Matrix{Float64,Float64} format
# -> the function returns also the center in Vector{Float64} format
function centering(points::Vector{Tuple{Float64,Float64}})
   n = length(points)
   center = zeros(2)
   for i = 1:n
      center[1] = center[1] + points[i][1]
      center[2] = center[2] + points[i][2]
   end
   center[1] = center[1] / n
   center[2] = center[2] / n

   # constructing the matrix of centered coordinates
   centered = zeros(2,n)
   for i = 1:n
      centered[1,i] = points[i][1] - center[1]
      centered[2,i] = points[i][2] - center[2]
   end

   # ending and returning the center
   return center, centered
end

# performing the Procrustes alignment of two point configurations
# -> it returns the aligned coordinates of points1 in Matrix{Float64,Float64} format
function align_procrustes(points1::Vector{Tuple{Float64,Float64}},points2::Vector{Tuple{Float64,Float64}})
   n = length(points1)
   if (n != length(points2)) throw(ArgumentError("The two input point configurations differ in length")) end

   # centering the data
   center1, centered1 = centering(points1)
   center2, centered2 = centering(points2)

   # computing the optimal rotation
   U,S,V = svd(centered2*centered1')
   R = V*U'

   # coomputing the optimal scaling
   scale = sum(centered2 .* centered2) / sum((R*centered1) .* (R*centered1))
   R = scale*R

   # computing the optimal translation
   translation = center2 - R*center1

   # applying the transformation and returning
   matrix = zeros(2,n)
   for i = 1:n
      matrix[1,i] = points1[i][1]
      matrix[2,i] = points1[i][2]
   end
   return R*matrix .+ translation
end

#  TEST0
#- we generate a random point configuration
#- we compute the distance and orientation matrix
#- we reconstruct the original configuration by invoking the 'reconstruct' function
function test0(n::Int64)
   println("Stick model test0")
   if (n <= 0) throw(ArgumentError("Problem size must be strictly positive")) end
   println("Problem size : ",n)
   points = gen_coordinates(n)
   distances = distance_matrix(points)
   orientations = orientation_matrix(points,distances)
   solution = reconstruct(distances,orientations)
   sol_distances = distance_matrix(solution)
   sol_orientations = orientation_matrix(solution,sol_distances)
   error = distance_error(distances,sol_distances)
   @printf "Error on distances : %20.16f\n" error
   error = orientation_error(orientations,sol_orientations)
   @printf "Error on orientations : %20.16f\n" error
end

#  TEST1
#- we generate a random point configuration
#- we compute the distance matrix and we perturb its entries
#- we compute the orientation matrix and we pertub the angles they represent
#- we remove the distance / orientation information for some pairs of points
#- we attempt reconstructing the original configuration through the Stick Model
function test1(n::Int64,epsilon::Float64,threshold::Float64)
   println("Stick model test1")
   if (n <= 0) throw(ArgumentError("Problem size must be strictly positive")) end
   println("Problem size : ",n)
   if (epsilon <= 0.0) throw(ArgumentError("Introduced epsilon error must be strictly positive")) end
   println("Introduced error : ",epsilon)
   if (threshold <= 0.0) throw(ArgumentError("Threshold value must be strictly positive")) end
   if (threshold >  1.0) throw(ArgumentError("Threshold value should not be larger than 1")) end
   println("Density threshold : ",threshold)

   # generating the point configuration and constructing problem instance
   points = gen_coordinates(n)
   distances = distance_matrix(points)
   perturb!(distances,epsilon)
   orientations = orientation_matrix(points,distances)
   perturb!(orientations,epsilon)
   sparsify!(distances,orientations,threshold)

   # generating and solving the Stick Model
   println("Solving the Stick Model ...")
   seconds = time()
   error, solution = stick_model(distances,orientations)
   seconds = time() - seconds
   @printf "Max error on stick coordinates : %20.16f\n" error
   sol_distances = distance_matrix(solution)
   error = distance_error(distances,sol_distances)
   @printf "Error on distances : %20.16f\n" error
   sol_orientations = orientation_matrix(solution,sol_distances)
   error = orientation_error(orientations,sol_orientations)
   @printf "Error on orientations : %20.16f\n" error
   @printf "Time : %20.16f\n" seconds

   # visualization
   original = zeros(2,n)
   for i = 1:n
      original[1,i] = points[i][1]
      original[2,i] = points[i][2]
   end
   aligned = align_procrustes(solution,points)
   plt = plot(original[1,:],original[2,:],seriestype=:scatter,label="Original Points",title="Original vs Reconstructed Points",xlabel="X",ylabel = "y",color="blue",markersize=5)
   plot!(plt,aligned[1,:],aligned[2,:],seriestype=:scatter,label="Reconstructed Points",color="red",markersize=5)  # solutions
   for i = 1:size(original,2)
      plot!(plt,[original[1,i],aligned[1,i]],[original[2,i],aligned[2,i]],color="gray",linewidth=0.5,label="")  # lines
   end
   display(plt)
end

#  TEST2
#- we generate a random point configuration
#- we compute the distance matrix and we perturb its entries
#- we randomly generate the orientation matrix
#- we remove the distance / orientation information for some pairs of points
#- we attempt reconstructing the original configuration through the Stick Heuristic
function test2(n::Int64,epsilon::Float64,threshold::Float64)
   println("Stick model test2")
   if (n <= 0) throw(ArgumentError("Problem size must be strictly positive")) end
   println("Problem size : ",n)
   if (epsilon <= 0.0) throw(ArgumentError("Introduced epsilon error must be strictly positive")) end
   println("Introduced error : ",epsilon)
   if (threshold <= 0.0) throw(ArgumentError("Threshold value must be strictly positive")) end
   if (threshold >  1.0) throw(ArgumentError("Threshold value should not be larger than 1")) end
   println("Density threshold : ",threshold)

   # generating the point configuration and constructing problem instance
   points = gen_coordinates(n)
   distances = distance_matrix(points)
   perturb!(distances,epsilon)
   orientations = random_orientations(n)
   sparsify!(distances,orientations,threshold)

   # solving the generated instance by the Stick Heuristic
   println("Running the Stick Heuristic ...")
   seconds = time()
   error, solution = stick_heuristic(2000,distances,0.01)
   seconds = time() - seconds
   @printf "Max error on stick coordinates : %20.16f\n" error
   sol_distances = distance_matrix(solution)
   error = distance_error(distances,sol_distances)
   @printf "Error on distances : %20.16f\n" error
   @printf "Time : %20.16f\n" seconds

   # visualization
   original = zeros(2,n)
   for i = 1:n
      original[1,i] = points[i][1]
      original[2,i] = points[i][2]
   end
   aligned = align_procrustes(solution,points)
   plt = plot(original[1,:],original[2,:],seriestype=:scatter,label="Original Points",title="Original vs Reconstructed Points",xlabel="X",ylabel = "y",color="blue",markersize=5)
   plot!(plt,aligned[1,:],aligned[2,:],seriestype=:scatter,label="Reconstructed Points",color="red",markersize=5)  # solutions
   for i = 1:size(original,2)
      plot!(plt,[original[1,i],aligned[1,i]],[original[2,i],aligned[2,i]],color="gray",linewidth=0.5,label="")  # lines
   end
   display(plt)
end

