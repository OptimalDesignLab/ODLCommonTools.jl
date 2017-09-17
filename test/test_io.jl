
facts("----- Testing IO -----") do

  fname = "test.dat"
  q = rand(Complex128, 2, 3, 4)
  write_binary(fname, q)

  q2 = zeros(q)
  read_binary!(fname, q2)

  @fact vecnorm(q - q2) --> roughly(0.0, atol=1e-14)

  # test error conditions
  q = zeros(Complex128, 2, 3)
  @fact_throws read_binary!(fname, q)


  q = zeros(Complex, 2, 2)
  @fact_throws write_binary(fname, q)
  @fact_throws read_binary!(fname, q)

  opts = Dict{Any, Any}()
  mesh = ODLCommonTools.ConcreteDGMesh{Float64}(opts)
  sbp = ODLCommonTools.ConcreteSBP()
  eqn = ODLCommonTools.ConcreteSolutionData{Float64, Float64}(mesh, sbp, opts)
end
