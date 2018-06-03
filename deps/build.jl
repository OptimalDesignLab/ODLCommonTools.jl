# Install dependencies

# install PkgFix if not present
if !isdir(joinpath(Pkg.dir(), "PkgFix"))
  Pkg.clone("https://github.com/OptimalDesignLab/PkgFix.jl.git")
end
Pkg.checkout("PkgFix", "upgrade_0.6")


using PkgFix  # from now on, use PkgFix instead of Pkg for everything



global const ARRAYVIEWS_URL = "https://github.com/JaredCrean2/ArrayViews.jl.git"
global const ARRAYVIEWS_VER = "work"

pkg_dict = PkgFix.installed()

if !haskey(pkg_dict, "ArrayViews")
  PkgFix.add(ARRAYVIEWS_URL, branch_ish=ARRAYVIEWS_VER)
else
  PkgFix.checkout("ArrayViews", ARRAYVIEWS_VER)
end


