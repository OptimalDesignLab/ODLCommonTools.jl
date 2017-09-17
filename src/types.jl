# sample implementations of the Abstract Types 
# used for testing, not to be used for real code

type ConcreteDGMesh{Tmsh} <: AbstractDGMesh{Tmsh}
  # counts
  numVert::Int
  numEl::Int
  numNodes::Int
  numDof::Int
  numDofPerNode::Int
  numNodesPerElement::Int
  order::Int
  numNodesPerFace::Int

  # parallel counts
  npeers::Int
  numGlobalEl::Int
  numSharedEl::Int
  peer_face_counts::Array{Int, 1}
  local_element_counts::Array{Int, 1}
  remote_element_counts::Array{Int, 1}

  # MPI info
  comm::Int  # this should really be an MPI.Comm
  myrank::Int
  commsize::Int
  peer_parts::Array{Int, 1}

  # Discretization type
  isDG::Bool
  isInterpolated::Bool

  # Mesh data
  coords::Array{Tmsh, 3}
  dxidx::Array{Tmsh, 4}
  jac::Array{Tmsh, 2}

  # interpolated data
  coords_bndry::Array{Tmsh, 3}
  dxidx_bndry::Array{Tmsh, 4}
  jac_bndry::Array{Tmsh, 2}
  dxidx_face::Array{Tmsh, 4}
  jac_face::Array{Tmsh, 2}

  # parallel data
  coords_sharedface::Array{Array{Tmsh, 3}, 1}
  dxidx_sharedface::Array{Array{Tmsh, 4}, 1}
  jac_sharedface::Array{Array{Tmsh, 2}, 1}

  # boundary condition data
  numBC::Int
  numBoundaryEdges::Int
  bndryfaces::Array{Boundary, 1}
  bndry_offsets::Array{Int, 1}
  bndry_funcs::Array{BCType, 1}

  # interior edge data
  numInterfaces::Int
  interfaces::Array{Interface, 1}

  # dof number data
  dofs::Array{Int, 3}
  dof_offset::Int
  sparsity_bnds::Array{Int, 2}
  sparsity_nodebnds::Array{Int, 2}

  # mesh coloring data
  numColors::Int
  maxColors::Int
  color_masks::Array{BitArray{1}, 1}
  shared_element_colormasks::Array{Array{BitArray{1}, 1}, 1}
  pertNeighborEls::Array{Int, 2}

  # parallel bookkeeping
  bndries_local::Array{Array{Boundary, 1}, 1}
  bndries_remote::Array{Array{Boundary, 1}, 1}
  shared_interfaces::Array{Array{Interface, 1}, 1}
  shared_element_offsets::Array{Int, 1}
  local_element_lists::Array{Array{Int, 1}, 1}

  function ConcreteDGMesh(opts)

    mesh = new()  # incomplete initialization

    # make up values
    mesh.numVert = 5
    mesh.numEl = 7
    mesh.numNodesPerElement = 4
    mesh.numNodes = mesh.numEl*mesh.numNodesPerElement
    mesh.numDofPerNode = 5
    mesh.numDof = mesh.numNodes*mesh.numDofPerNode
    mesh.order = 2
    mesh.numNodesPerFace = 3

    mesh.npeers = 2
    mesh.numGlobalEl = 42
    mesh.numSharedEl = 4
    mesh.peer_face_counts = [2, 2]
    mesh.local_element_counts = [2, 2]
    mesh.remote_element_counts = [1, 2]

    mesh.comm = 666
    mesh.myrank = 0
    mesh.commsize = 4
    mesh.peer_parts = [1, 3]

    mesh.isDG = true
    mesh.isInterpolated = true

    mesh.coords = rand(2, mesh.numNodesPerElement, mesh.numEl)
    mesh.dxidx = rand(2, 2, mesh.numNodesPerElement, mesh.numEl)
    mesh.jac = rand(mesh.numNodesPerElement, mesh.numEl)

    mesh.numBoundaryEdges = 2
    mesh.coords_bndry = rand(2, mesh.numNodesPerFace, mesh.numBoundaryEdges)
    mesh.dxidx_bndry = rand(2, 2, mesh.numNodesPerFace, mesh.numBoundaryEdges)
    mesh.jac_bndry = rand(mesh.numNodesPerFace, mesh.numBoundaryEdges)

    mesh.numInterfaces = 2
    mesh.dxidx_face = rand(2, 2, mesh.numNodesPerFace, mesh.numInterfaces)
    mesh.jac_face = rand(mesh.numNodesPerFace, mesh.numInterfaces)

    mesh.numBC = 2
    mesh.bndryfaces = [Boundary(2, 1), Boundary(3, 2)]
    mesh.bndry_offsets = [1, 2, 3]
    mesh.bndry_funcs = Array(BCType, mesh.numBC)

    mesh.interfaces = [Interface(2, 3, 2, 3, 0), Interface(4, 7, 1, 3, 2)]

    mesh.dofs = Array(Int, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    println("numDof = ", mesh.numDof)
    println("size(mesh.dofs) = ", size(mesh.dofs))

    for i=1:mesh.numDof
      mesh.dofs[i] = i
    end

    mesh.dof_offset = 0
    mesh.sparsity_bnds = Array(Int, 2, mesh.numDof)
    # make the matrix dense
    for i=1:mesh.numDof
      mesh.sparsity_bnds[1, i] = 1
      mesh.sparsity_bnds[2, i] = mesh.numDof
    end

    mesh.sparsity_nodebnds = Array(Int, 2, mesh.numNodes)
    for i=1:mesh.numNodes
      mesh.sparsity_nodebnds[1, i] = 1
      mesh.sparsity_nodebnds[2, i] = mesh.numNodes
    end

    mesh.numColors = 4  # like that's really possible...
    mesh.maxColors = 5
    mesh.color_masks = Array(BitArray{1}, mesh.numColors)
    mesh.shared_element_colormasks = Array(Array{BitArray{1}, 1}, mesh.npeers)
    mesh.pertNeighborEls = Array(Int, 4, mesh.numEl)

    mesh.bndries_local = Array(Array{Boundary, 1}, mesh.npeers)
    mesh.bndries_remote = Array(Array{Boundary, 1}, mesh.npeers)
    mesh.shared_interfaces = Array(Array{Interface, 1}, mesh.npeers)
    mesh.shared_element_offsets = Array(Int, 2)
    mesh.local_element_lists = Array(Array{Int, 1}, mesh.npeers)

    return mesh
  end

end



type ConcreteParamType{Tdim}
  t::Float64
  order::Int
  time  # this really should have a concrete type

  function ConcreteParamType(mesh, opts)

    t = 0.0
    order = mesh.order
    time = 0

    return new(t, order, time)
  end
end

type ConcreteSolutionData{Tsol, Tres} <: AbstractSolutionData{Tsol, Tres}
  q::Array{Tsol, 3}
  q_vec::Array{Tsol, 1}
  # AbstractSharedFaceData should be a concrete type
  shared_data::Array{AbstractSharedFaceData{Tsol}, 1}
  res::Array{Tres, 3}
  res_vec::Array{Tres, 1}
  M::Array{Float64, 1}
  Minv::Array{Float64, 1}
  disassembleSolution::Function
  assembleSolution::Function
  multiplyA0inv::Function
  majorIterationCallback::Function
  params::ConcreteParamType{2}

  function ConcreteSolutionData(mesh, sbp, opts)

    numDofPerNode = mesh.numDofPerNode
    numNodesPerElement = mesh.numNodesPerElement
    numEl = mesh.numEl
    numNodesPerFace = mesh.numNodesPerFace
    numDof = mesh.numDof

    q = Array(Tsol, numDofPerNode, numNodesPerElement, numEl)
    q_vec = Array(Tsol, numDof)
    # the elements of this array should be populated as well
    shared_data = Array(AbstractSharedFaceData{Tsol}, mesh.npeers)
    res = Array(Tres, mesh.numDofPerNode, mesh.numNodesPerElement, mesh.numEl)
    res_vec = Array(Tres, mesh.numDof)
    M = ones(numDof)
    Minv = ones(numDof)
    
    # create empty functions
    disassembleSolution = (mesh, sbp, eqn, opts) -> begin end
    assembleSolution = (mesh, sbp, eqn, opts) -> begin end
    multiplyA0inv = (mesh, sbp, eqn, opts, res_arr) -> begin end
    majorIterationCallback = (itr, mesh, sp, eqn, opts) -> begin end

    params = ConcreteParamType{2}(mesh, opts)

    return new(q, q_vec, shared_data, res, res_vec, M, Minv,
               disassembleSolution, assembleSolution, multiplyA0inv,
               majorIterationCallback, params)
  end
end

type ConcreteSBP
  Q::Array{Float64, 2}

  function ConcreteSBP()
    Q = rand(3, 3)
    return new(Q)
  end
end
