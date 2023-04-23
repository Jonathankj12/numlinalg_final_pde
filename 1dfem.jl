using LinearAlgebra
using Plots
using Integrals
using QuadGK
using ForwardDiff


l=3.0 # domain length
m=30 # number of nodes m-1 is number of elements
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

n_hat(i,m,x)=hat_func(i,x,m)

nodes[2]


function sin_basis(k::Int,mesh,x)
    return sin((k-1)*pi*x/mesh[end])
end

function poly_basis(k::Int,mesh,x)
    return (x^k)*(x-mesh[end])
end


poly_basis(2,nodes,0)

sin_basis(2,nodes,3)

sip(x) = sin_basis(2,nodes,x)

plb(x)=poly_basis(2,nodes,x)

nh(x)=n_hat(2,nodes,x)

sip(0)


ForwardDiff.derivative(sip,3.0)

ForwardDiff.derivative(plb,3)

ForwardDiff.derivative(nh,nodes[2]*3/2)

function bi_form_arb(i,j,mesh,func_gen)
    f_basis_i(x) = func_gen(i,mesh,x)
    f_basis_j(x) = func_gen(j,mesh,x)
    Df_i(x)= ForwardDiff.derivative(f_basis_i,x)
    Df_j(x)= ForwardDiff.derivative(f_basis_j,x)
    integrand(x,p) = Df_i(x)*Df_j(x)
    prob=IntegralProblem(integrand,mesh[1],mesh[end])
    solve(prob,QuadGKJL()).u
end


bi_form_arb(4,7,nodes,sin_basis)

function lin_form_arb(i,mesh,func_gen,force)
    f_basis(x) = func_gen(i,mesh,x)
    integrand(x,p) = f_basis(x)*force(x)
    prob=IntegralProblem(integrand,mesh[1],mesh[end])
    solve(prob,QuadGKJL()).u
end






#rhs of the original poission equation, constant for simplicity
function init_force(x)
    return 3
end

lin_form_arb(3,nodes,sin_basis,init_force)

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
#A=[bilinear_form_pwlinmat(r,c,hs) for r =2:(m-1), c=2:(m-1)]

#calculate rhs of system elements acutally integrates the basis funcitons so could use different shape function
function lin_form_pwlinvec(i,node)
    integrand(x,p) = hat_func(i,x,node)*init_force(x)
    prob=IntegralProblem(integrand,node[i-1],node[i+1])
    res=solve(prob,QuadGKJL()).u
    return res
end

A=[bi_form_arb(r,c,nodes,poly_basis) for r =2:(m-1), c=2:(m-1)]
#construct rhs
#F= [lin_form_pwlinvec(i,nodes) for i= 2:(m-1)]
F= [lin_form_arb(i,nodes,poly_basis,init_force) for i= 2:(m-1)]
#solve system
coeffs=A\F

#construct finite element approximation using the shape funciton and the calculated coeffs
function U_h(coe,node,bfunc,x)
    basis= [bfunc(i,node,x) for i= 2:(m-1)]
    return u=dot(coe,basis)
end






#Plot the calculated and true solutions
xp=collect(range(0,l,300))
uh=U_h.(Ref(coeffs),Ref(nodes),poly_basis,xp)
f=3
utrue(x)=-((f/2)*x^2)+(((f*l)/2)*x)
ut=utrue.(xp)
plot(xp,[uh,ut])












