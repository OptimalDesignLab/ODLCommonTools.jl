# test ro_view
#using ODLCommonTools
#using FactCheck

function sum_cols{T}(A::AbstractMatrix{T})

  val = zero(T)
  for i=1:size(A, 2)
    tmp = ro_sview(A, :, i)
    val += sum(tmp)
  end

  return val
end

function sumit{T}(A::AbstractArray{T})

  val = zero(T)
  for i=1:length(A)
    val += A[i]
  end

  return val
end

facts("----- Testing ROView -----") do

  a = rand(100, 1000)

  b = ROView(a)
  @fact size(b) --> size(a)
  @fact length(b) --> length(a)

  @fact norm(a - b) --> roughly(0.0)

  # check that allocating a ROView does not allocate 
  if !ODLCommonTools.safe_views
    @time sum_cols(a)
    bytes_alloc = @allocated sum_cols(a)

    @fact bytes_alloc --> less_than(100)
  end

  # check that indexing a ROView does not allocate
  if Base.JLOptions().can_inline == true
    sumit(b)
    bytes_alloc = @allocated sumit(b)
    @fact bytes_alloc --> less_than(100)
  end

end

