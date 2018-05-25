using ODLCommonTools
using FactCheck
import ODLCommonTools.get_uninitialized_SolutionData

mutable struct AbstractSolutionData_ <: AbstractSolutionData
  q::Array{Float64,2}
  q_vec::Array{Float64,1}
  res::Array{Float64,2}
  res_vec::Array{Float64,1}
  shared_data::Array{Array{Float64, 3}, 1}  # this isn't exactly like real 
                                            # equation objects, but close enough
  unhandled_field::Float64
  function AbstractSolutionData_(eqn::AbstractSolutionData)
    return new()
  end

  
  function AbstractSolutionData_(
    q::Array{Float64,2}, 
    q_vec::Array{Float64,1},
    res::Array{Float64,2},
    res_vec::Array{Float64,1},
    shared_data::Array{Array{Float64,3}, 1},
    unhandled_field::Float64,
    )
    return new(q, q_vec, res, res_vec, shared_data, unhandled_field)
  end
  
#   function AbstractSolutionData_(
#     eqn = new()
#     eqn.q = q
#     eqn.q_vec = q_vec
#     eqn.res = res
#     eqn.res_vec = res_vec
#     eqn.q_face_send = q_face_send
#     eqn.q_face_recv = q_face_recv
#     eqn.unhandled_field = unhandled_field
#     return eqn
#     return new(q, q_vec, res, res_vec, q_face_send, q_face_recv, unhandled_field)
#   end

end

function get_uninitialized_SolutionData(eqn::AbstractSolutionData_)

  return AbstractSolutionData_(eqn)

end

function make_vector(a::AbstractArray{T}) where T

  v =zeros(T, length(a))
  for i=1:length(a)
    v[i] = a[i]
  end

  return v
end


facts("--- Testing copy_array_recursive ---") do

  # test single level
  arr = rand(5)
  arr_dest = zeros(5)
  ODLCommonTools.copy_array_recursive!(arr_dest, arr)
  @fact norm(arr_dest - arr) --> roughly(0.0, atol=1e-14)

  # test two levels
  arr = Array(Array{Float64, 2}, 2)
  arr_dest = Array(Array{Float64, 2}, 2)
  for i=1:2
    arr[i] = rand(i+2, i+2)
    arr_dest[i] = zeros(i+2, i+2)
  end

  ODLCommonTools.copy_array_recursive!(arr_dest, arr)
  for i=1:2
    @fact norm(arr_dest[i] - arr[i]) --> roughly(0.0, atol=1e-14)
  end

  # test three levels
  arr = Array(Array{Array{Float64, 1}, 1}, 2)
  arr_dest = Array(Array{Array{Float64, 1}, 1}, 2)
  for i=1:2
    arr[i] = Array(Array{Float64, 1}, 2)
    arr_dest[i] = Array(Array{Float64, 1}, 2)
    for j=1:2
      arr[i][j] = rand(3)
      arr_dest[i][j] = zeros(3)
    end
  end

  ODLCommonTools.copy_array_recursive!(arr_dest, arr)
  for i=1:2
    for j=1:2
      @fact norm(arr_dest[i][j] - arr[i][j]) --> roughly(0.0, atol=1e-14)
    end
  end

end

facts("--- Testing copy functions ---") do

  num_inner_arrays = 3

  array_16 = [[1.0, 2, 3] [4.0, 5, 6]]
  array_16_double = [[1.0, 2, 3] [4.0, 5, 6]]*2.0
  array_163_inner = reshape([1.0 2 3 4 5 6 7 8], 2, 2, 2)
  array_163 = Array(Array{Float64,3}, num_inner_arrays)
#  array_163_double = Array(Array{Float64,3}, num_inner_arrays)
  for i = 1:num_inner_arrays
    array_163[i] = array_163_inner
#    array_163_double[i] = array_163_inner*2.0
  end
  
#   println(array_163)

  test_data = AbstractSolutionData_(
    array_16,
    make_vector(array_16),
    array_16_double,
    make_vector(array_16),
    array_163,
 #   array_163_double,
    7.0,
  )

  mat1 = rand(3,2)
  mat2 = rand(3,2)

  mat_2223_inner = rand(2,2,2)
  mat_2223 = Array(Array{Float64,3}, num_inner_arrays)
#  mat_2223_double = Array(Array{Float64,3}, num_inner_arrays)
  for i = 1:num_inner_arrays
    mat_2223[i] = mat_2223_inner
#    mat_2223_double[i] = mat_2223_inner*2.0
  end
  float_field = rand(1)[1]

  test_data_rand = AbstractSolutionData_(
    mat1,
    make_vector(mat1),
    mat2,
    make_vector(mat2),
    mat_2223,
#    mat_2223_double,
    float_field
  )

#   println(test_data.q)
#   println(test_data_rand.q)
#   println(test_data.res)
#   println(test_data_rand.res)
#   println("---now copying---")

  #---------------
  # testing copy!(eqn_dest::AbstractSolutionData, eqn_src::AbstractSolutionData)
  copy!(test_data, test_data_rand)

#   println(test_data.q)
#   println(test_data_rand.q)
#   println(test_data.res)
#   println(test_data_rand.res)

  @fact test_data.q --> roughly(test_data_rand.q, atol=1e-14)
  @fact test_data.q_vec --> roughly(test_data_rand.q_vec, atol=1e-14)
  @fact test_data.res --> roughly(test_data_rand.res, atol=1e-14)
  @fact test_data.res_vec --> roughly(test_data_rand.res_vec, atol=1e-14)
  for i = 1:num_inner_arrays
    @fact test_data.shared_data[i] --> roughly(test_data_rand.shared_data[i], atol=1e-14)
#    @fact test_data.q_face_send[i] --> roughly(test_data_rand.q_face_send[i], atol=1e-14)
#    @fact test_data.q_face_recv[i] --> roughly(test_data_rand.q_face_recv[i], atol=1e-14)
  end

  #---------------
  # testing that eqn_dest = copy(eqn_src::AbstractSolutionData) throws an error
  @fact_throws test_data_dest = copy(test_data)

  #---------------
  # testing copyForMultistage!
  # use multistage copy eqn1 to eqn2, then write to eqn1, check that eqn2 fields are changed or unchanged as appropriate
  fill!(test_data.q, 0.0)
  fill!(test_data.q_vec, 0.0)
  fill!(test_data.res, 0.0)
  fill!(test_data.res_vec, 0.0)
  for i=1:length(test_data.shared_data)
    fill!(test_data.shared_data[i], 0.0)
  end
  test_data.unhandled_field = 17.0
  copyForMultistage!(test_data, test_data_rand)

#  test_data_copy = copyForMultistage!(test_data)
  @fact test_data.q --> roughly(test_data_rand.q, atol=1e-14)
  @fact test_data.q_vec --> roughly(test_data_rand.q_vec, atol=1e-14)
  @fact test_data.res --> roughly(test_data_rand.res, atol=1e-14)
  @fact test_data.res_vec --> roughly(test_data_rand.res_vec, atol=1e-14)
  for i = 1:num_inner_arrays
    @fact test_data.shared_data[i] --> roughly(test_data_rand.shared_data[i], atol=1e-14)
#    @fact test_data.q_face_send[i] --> roughly(test_data_rand.q_face_send[i], atol=1e-14)
#    @fact test_data.q_face_recv[i] --> roughly(test_data_rand.q_face_recv[i], atol=1e-14)
  end
#   println(test_data.unhandled_field)
#   println(test_data_rand.unhandled_field)
  @fact test_data.unhandled_field --> not(roughly(test_data_rand.unhandled_field, atol=1e-14))

  # verfy no aliasing
  @fact pointer(test_data.q) --> not( pointer(test_data_rand.q) )
  @fact pointer(test_data.shared_data) --> not( pointer(test_data_rand.q) )
  for i=1:length(test_data.shared_data)
    @fact pointer(test_data.shared_data[i]) --> not( pointer(test_data_rand.shared_data[i]))
  end
#   println(test_data)
#   println(test_data_copy)

end


