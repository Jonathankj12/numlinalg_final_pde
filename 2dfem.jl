using LinearAlgebra
using Plots
using Integrals
import Gmsh: gmsh

function create_fem_mesh()
    gmsh.initialize()
    gmsh.model.add("2d_fem_model")
    lc = 1e-2
    gmsh.model.geo.addPoint(0,0,0,lc,1)
    gmsh.model.geo.addPoint(.1,0,0,lc,2)
    gmsh.model.geo.addPoint(.1,.3,0,lc,3)
    p4 = gmsh.model.geo.addPoint(0,.3,0,lc)
    gmsh.model.geo.addLine(1,2,1)
    gmsh.model.geo.addLine(3,2,2)
    gmsh.model.geo.addLine(3,p4,3)
    gmsh.model.geo.addLine(4,1,p4)
    gmsh.model.geo.addCurveLoop([4,1,-2,3],1)
    gmsh.model.geo.addPlaneSurface([1],1)
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(0,[1,2],1)
    gmsh.model.addPhysicalGroup(1,[1,2],2)
    gmsh.model.addPhysicalGroup(2,[1],6)
    gmsh.model.setPhysicalName(2,6,"My surface")
    mesh=gmsh.model.mesh.generate(2)
    #gmsh.write("t1.msh")
    return mesh
end

create_fem_mesh()


gmsh.isInitialized()

mutable struct domain_mesh
end

function clean_up_fem()
    if (gmsh.isInitialized() == 1)
        gmsh.finalize()
    end
end

clean_up_fem()

