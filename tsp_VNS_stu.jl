# -------------------------------------------------------------------------------------------
# X. Gandibleux - Metaheuristiques : TSP - Septembre 2019
# Base de depart « etudiant » pour exercice "VNS applique au TSP"

@static if VERSION < v"1.0-"
    error("NOT COMPLIANT WITH JULIA < v1.0.0")
end

using Random, Printf

# -------------------------------------------------------------------------------------------
# construit un distancier symetrique aux valeurs de distances uniformément aléatoires dans 1..vSup
function tsp_elaboreDistancierAleatoire(
           n::Int64,    # in : nombre de villes
           vSup::Int64  # in : plus grande distance entre deux villes
         )
    distancier = rand(1:vSup,n,n)
    for i=1:n
        for j=1:n
            distancier[i,j] = distancier[j,i]
        end
        distancier[i,i] = 0
    end
    return distancier # out : matrice carree symetrique de distances entieres
end

# -------------------------------------------------------------------------------------------
# calcule la longueur d’un circuit hamiltonien correspondant a une solution TSP admissible
function tsp_calculeLongueurCircuit(
           d::Matrix{Int64},  # in : distancier symetrique de ville a ville
           x::Array{Int64},   # in : solution admissible représentée par codage de permutation des villes
         )
    n::Int64 = length(x) # Nombre de villes
    longueur::Int64 = d[ x[n] , x[1] ]
    for i = 1:n-1
        longueur += d[ x[i] , x[i+1] ]
    end
    return longueur # out : longueur du circuit correspondant a la solution admissible
end

# -------------------------------------------------------------------------------------------
# construit une solution aléatoire TSP admissible représentée dans un codage de permutation
function tsp_elaboreSolutionAleatoire(
           n::Int64  # in : nombre de villes
         )
    return shuffle(collect(1:n))  # out : une permutation aleatoire des villes
end

# -------------------------------------------------------------------------------------------
# echange deux valeurs dans un vecteur sur des positions tirees aleatoirement
function swap!(
           x::Array{Int64}   # inout : vecteur d’entiers (codage binaire ou codage de permutation)
         )
    # swap entre informations aux positions i1 et i2
    i1 = rand(1:length(x))
    i2 = rand(1:length(x))
    x[i1], x[i2] = x[i2], x[i1]
end

# -------------------------------------------------------------------------------------------
function tsp_AppliqueLocalSearch2opt(
           d::Matrix{Int64}, # in : distancier symetrique de ville a ville
           x::Array{Int64},  # in : solution TSP admissible représentée pas un codage de permutation des villes
           z::Int64          # in : valeur de la solution TSP admissible
         )
    # A mettre en place ...
    return xSol, zSol    # out : meilleure solution TSP admissible trouvee et sa valeur
end

# -------------------------------------------------------------------------------------------
function tsp_solveVNS(
           d::Matrix{Int64},     # in : distancier symetrique de ville a ville
           xBest::Array{Int64},  # in : solution TSP admissible représentée pas un codage de permutation des villes
           zBest::Int64          # in : valeur de la solution TSP admissible
         )
    # A mettre en place ...
    return xBest, zBest, pBest    # out : meilleure solution TSP admissible trouvee, sa valeur, la courbe
end

# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------

function tsp_VNS(
           n::Int64    = 20,  # in : nombre de villes
           vSup::Int64 = 25,  # in : plus grande distance entre deux villes
           nRun::Int64 = 1,  # in : nombre de run VNS a realiser
         )

    d = tsp_elaboreDistancierAleatoire(n,vSup) # distancier aleatoire

    for run =1:nRun
        xInit = tsp_elaboreSolutionAleatoire(n)
        zInit = tsp_calculeLongueurCircuit(d,xInit)

        println("Init    : z(xInit) = ", zInit)
        #@show xInit
        println("  ")

        xBest, zBest, pBest = tsp_solveVNS(d,xInit,zInit)
        println("VNS     : z(xVNS)  = ", zBest)
        #@show xBest
        println("  ")
    end
end
