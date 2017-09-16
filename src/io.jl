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
function write_binary{T}(fname::AbstractString, arr::ContiguousArrays{T})

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
  Reads files written by [`write_binary1}(@ref)

  **Inputs**

   * fname: file name, including extension.  Can be relative of absolute path

  **Inputs/Outputs*

   * arr: array to be populated with the contents of the file.  The array
          must be the same size as the file.
"""
function read_binary!{T}(fname::AbstractString, arr::ContiguousArrays{T})

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
  This function writes the solution vector (eqn.qvec) to a file.
  The path to the file must already exist.  If the file already exists, it
  will be overwritten.

  **Inputs**

   * mesh: an AbstractMesh
   * sbp: an AbstractSBP
   * eqn: an AbstractSolutionData
   * opts: options dictionary
   * fname: file name, without extension.  Can be absolute or relative path.
            In parallel, the file name should *not* contain the MPI rank.

"""
function writeSolutionFiles(mesh::AbstractMesh, sbp, eqn::AbstractSolutionData,
                            opts, fname::AbstractString)


  # create parallel file name
  
  # add extension

  # call write_binary
 
  return nothing
end

#TODO: create a CheckPointer type and functions that manage it

