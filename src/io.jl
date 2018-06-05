# IO for the AbstractSolutionData (and maybe other?) objects

"""
  This function writes an array to a file an an efficient manner.  Ovewrites
  existing file with same name.

  **Inputs**

   * fname: the file name, including extension.  Can be a relative or absolute
            path
   * arr: the array to write.  Its memory representation must be contiguous
         

  Implementation:
    Currently this function calls julias `write` function, but in the future
    we may switch to a more elaborate file format.  Any machinery that is
    built on top of this should not need to know the file format, hence only
    this function will need to be modified (and the associated reader function)
"""
function write_binary(fname::AbstractString, arr::ContiguousArrays{T}) where T

  if !isbits(T)
    throw(ErrorException("cannot write array of non-isbits elements to binary file"))
  end

  f = open(fname, "w")
  n = write(f, arr)
  close(f)

  @assert n == sizeof(arr)
  
  return n
end

"""
  Reads files written by [`write_binary`](@ref).

  **Inputs**

   * fname: file name, including extension.  Can be relative of absolute path

  **Inputs/Outputs*

   * arr: array to be populated with the contents of the file.  The array
          must be the same size as the file.


  Implementation notes:

    The load will succeed as long as the array is the same size as the file,
    even if the datatype is different.

"""
function read_binary!(fname::AbstractString, arr::ContiguousArrays{T}) where T

  if !isbits(T)
    throw(ErrorException("cannot read array of non-isbits elements to binary file"))
  end

  f = open(fname, "r")

  # check the file is the right size
  fsize = filesize(f)
  if fsize != sizeof(arr)
    close(f)  # make sure to close the file if the exception will be thrown
    @assert fsize == sizeof(arr)
  end

  read!(f, arr)
  close(f)

end

"""
  This function writes the solution vector (eqn.q_vec) to a file.
  The path to the file must already exist.  If the file already exists, it
  will be overwritten.

  This function should be used as though it is collective on mesh.comm,
  although the current implementation may or may not require it.

  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an AbstractSBP
   * eqn: an AbstractSolutionData
   * opts: options dictionary
   * fname: file name, without extension.  Can be absolute or relative path.
            In parallel, the file name should *not* contain the MPI rank
            (this function will add it internally)
"""
function writeSolutionFiles(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData,
                            opts, fname::AbstractString)


  # create parallel file name
  fname = get_parallel_fname(fname, mesh.myrank)
  # add extension
  fname = string(fname, ".dat")

  # call write_binary
  write_binary(fname, eqn.q_vec)
 
  return nothing
end

"""
  This function reads the files written by [`writeSolutionFiles`](@ref)
 
  This function should be used as though it is collective on mesh.comm,
  although the current implementation may or may not require it.

 
  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an AbstractSBP
   * opts: options dictionary
   * fname: file name, without extension.  Can be absolute or relative path.
            In parallel, the file name should *not* contain the MPI rank

  **Inputs/Outputs**

   * eqn: eqn.q_vec is overwritten with the data loaded from the file.
          The vector must be the same length and have the same element type
          as the vector that was saved to the file.  No error is given if
          these conditions are not met
   
"""
function readSolutionFiles(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData,
                            opts, fname::AbstractString)


  # create parallel file name
  fname = get_parallel_fname(fname, mesh.myrank)
  
  # add extension
  fname = string(fname, ".dat")

  # call write_binary
  read_binary!(fname, eqn.q_vec)

  return nothing
end
