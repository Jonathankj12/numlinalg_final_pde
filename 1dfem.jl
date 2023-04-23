using LinearAlgebra
using Plots
using Integrals
using QuadGK


l=3.0 # domain length
m=10 # number of nodes m-1 is number of elements
nodes=collect(range(0,l,m)) #comment this out to test with non-uniform mesh

#Manually placed node points
#nodes=zeros(Float64,m)
#nodes[m]=l
#nodes[1]=0.0
#nodes[2]=l/2
#nodes[3]=3*l/4

#Nodes that half the distance between them each time
#nodes[1]=l/2
#for k in 2:m-1
#    nodes[k]=nodes[k-1]+(l/(2^k))
#end
#nodes[m]=l



#Element lengths
hs=zeros(Float64,m-1)
for i in 1:(m-1)
    hs[i]= nodes[i+1]-nodes[i]
end

#Element function hat
function hat_func(index::Int,x,mesh)
    if (index <= 1) || (index >= m)
        throw(DomainError(index,"basis functions are only defined on the interior of the domain"))
    end
    step1=mesh[index]-mesh[index-1]
    step2=mesh[index+1]-mesh[index]
    out = 0.0
    if (mesh[index-1] <= x) && (x < mesh[index])
        out=(x-mesh[index-1])/step1
    elseif ((mesh[index] <= x) && (x < mesh[index+1]))
        out=(mesh[index+1]-x)/step2
    else
        out=0
    end
    return out
end

#rhs of the original poission equation, constant for simplicity
function init_force(x)
    return 3
end

#Calculates the stiffness matrix elements, currently ues knowlege of hat function derivative, for simplicity.
function bilinear_form_pwlinmat(i,j,elems)
    res=0.0
    if i == j
        res=(1/elems[i-1]) + (1/elems[i])
    elseif (i-j) == 1 
        res=(-1/elems[i-1])
    elseif (j-i) == 1
        res = (-1/elems[j-1])
    else
        res = 0.0
    end
    return res
end

#construnct matrix
A=[bilinear_form_pwlinmat(r,c,hs) for r =2:(m-1), c=2:(m-1)]

#calculate rhs of system elements acutally integrates the basis funcitons so could use different shape function
function lin_form_pwlinvec(i,node)
    integrand(x,p) = hat_func(i,x,node)*init_force(x)
    prob=IntegralProblem(integrand,node[i-1],node[i+1])
    res=solve(prob,QuadGKJL()).u
    return res
end

#construct rhs
F= [lin_form_pwlinvec(i,nodes) for i= 2:(m-1)]
#solve system
coeffs=A\F

#construct finite element approximation using the shape funciton and the calculated coeffs
function U_h(coe,x,node)
    basis= [hat_func(i,x,node) for i= 2:(m-1)]
    return u=dot(coe,basis)
end

#Plot the calculated and true solutions
xp=collect(range(0,l,30))
uh=U_h.(Ref(coeffs),xp,Ref(nodes))
f=3
utrue(x)=-((f/2)*x^2)+(((f*l)/2)*x)
ut=utrue.(xp)
plot(xp,[uh,ut])












