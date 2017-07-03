# types and functions describing the topology of the reference element

typealias TopoIdxType Int

"""
### ODLCommonTools.ElementTopology

  This type describes the topology of the reference element.  For now, it only
  contains face information, but eventually will needed more info.


  The way to think about element topology is that an element is defined by its
  vertices, numbered 1 through n.  All other parts of the element's topology
  (edges, faces, etc.) can be defined in terms of the vertices.  It is then
  possible to construct the relationship between the higher dimensional
  topological entities directly (ie. edges to faces, etc.), although this is
  not yet implemented.

  Terminology:
    numFaces: the number of faces an element has.  A triangle has 3 faces (the
              edges) and a tet has 4 faces
    numEdges: The number of edges 

  TODO: consider making the indices Cints or UInt8s?

  Fields:
    face_verts: Tdim x numFaces array holding the indices of the vertices of
                each face
    edge_verts: The 2 x numEdges array holding the indices of the vertices of
                each face.  Note that this array specifies the orientation of
                the edge, ie. the edge is directed from the  first vertex to
                the second vertex.
                In 2D this is the same as face_verts.

    face_edges: Tdim x numEdgesPerFace x numFaces array containing the
                indices of the vertices (first dimension) that define each edge
                (second dimension) of each face (3rd dimension)

  Constructor:
    face_verts: the indices of the vertices that define each face
    edge_verts: the indices of the vertices that define each edge.
                If there are not know, thie argument can be omitted, however
                the face_edges field will not be usable
    topo2:      Only needed in 3D when edge_verts are supplied.  Describes the
                topology of a face of the element
"""

immutable ElementTopology{Tdim}
  face_verts::Array{TopoIdxType, 2}
  edge_verts::Array{TopoIdxType, 2}
  face_edges::Array{TopoIdxType, 2}
  face_edges_flipped::Array{Bool, 2}

  function ElementTopology{I1 <: Integer, I2 <: Integer}(
                           face_verts::AbstractArray{I1, 2}, 
                           edge_verts::AbstractArray{I2, 2}=zeros(Int, 0,0);
                           topo2::ElementTopology=ElementTopology{1}() )

    # do sanity checks

    # check all vertices are within range
    for i=1:length(face_verts)
      @assert face_verts[i] > 0
      @assert face_verts[i] <= Tdim + 1
    end

    # check faces are distinct
    println("face_verts = \n", face_verts)
    for i=1:size(face_verts, 2)
      curr_face = sort(face_verts[:, i])
      for j=(i+1):size(face_verts, 2)
        @assert sort(face_verts[:, j]) != curr_face
      end
    end

    if Tdim == 2
      edge_verts = face_verts
      face_edges = Array(TopoIdxType, 0, 0)
      face_edges_flipped = Array(Bool, 0, 0)
    elseif Tdim == 3
      if size(edge_verts, 1) == 0  # edge verts not known
        face_edges = Array(TopoIdxType, 0, 0)
        face_edges_flipped = Array(Bool, 0, 0)
      else
        face_edges, face_edges_flipped = construct_face_edges(topo2, face_verts, edge_verts)
      end
    end

    return new(face_verts, edge_verts, face_edges, face_edges_flipped)

  end  # end function

  function ElementTopology()
    face_verts = zeros(TopoIdxType, 0, 0)
    edge_verts = zeros(face_verts)
    face_edges = zeros(face_verts)
    face_edges_flipped = Array(Bool, 0, 0)
    return new(face_verts, edge_verts, face_edges, face_edges_flipped)
  end
end

"""
  Default constructor that uses Pumi topology
"""
function ElementTopology3()
  face_verts = [1 1 2 1; 2 2 3 3; 3 4 4 4]
  return ElementTopology{3}(face_verts)
end

"""
  Constructs a dummy 2d ElementTopology
"""
function ElementTopology2()
  face_verts = [1 2 3; 2 3 1]
  return ElementTopology{2}(face_verts)
end


#------------------------------------------------------------------------------
# functions to assist the constructor

"""
  This function uses a brute force lookup approach to compute the
  face_edges array.  Throws an exception if the algorithm cannot construct
  the array.

  Inputs:
    topo2: an ElementTopology2 to related the vertices of each face to their
           edges
    face_verts:
    edge_verts:

  Outputs:
    
"""
function construct_face_edges(topo2::ElementTopology{2},
                              face_verts::AbstractMatrix,
                              edge_verts::AbstractMatrix)
  dim = size(face_verts, 1)
  numFacesPerElement = size(face_verts, 2)
  numEdgesPerElement = size(edge_verts, 2)
  numEdgesPerFace = size(topo2.face_verts, 2)

  # index of the edge of each face in edge_verts
  face_edges = zeros(TopoIdxType, numEdgesPerFace, numFacesPerElement)

  # is the edge of the face oriented the same as the edge of the tet
  face_edges_flipped = Array(Bool, numEdgesPerFace, numFacesPerElement)

  for i=1:numFacesPerElement
    for j=1:numEdgesPerFace  # loop over edges on this face
      v1 = face_verts[topo2.face_verts[1, j], i]
      v2 = face_verts[topo2.face_verts[2, j], i]

      edge_verts_k = [v1, v2]

      # search edge list 
      is_flipped = false
      matched_edgenum = 0
      for k=1:numEdgesPerElement
        matches_edge, is_flipped = compare_edge(edge_verts_k, edge_verts[:, k])
        if matches_edge
          matched_edgenum = k; break
        end
      end

      @assert matched_edgenum != 0  # we should have found the edge somewhere

      face_edges[j, i] = matched_edgenum
      face_edges_flipped[j, i] = is_flipped

    end
  end

  return face_edges, face_edges_flipped
end


function compare_edge(edges1::AbstractVector, edges2::AbstractVector)

  matches = false
  is_flipped = false  # value doesn't matter if matches is false
  if edges1[1] == edges2[1] && edges1[2] == edges2[2]
    matches = true
    is_flipped = false
  elseif edges1[1] == edges2[2] && edges1[2] == edges2[1]
    matches = true
    is_flipped = true
  else
    matches = false
    is_flipped = false  # arbitrary
  end

  return matches, is_flipped
end

