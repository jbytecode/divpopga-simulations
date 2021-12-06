# CEC2017 SINGLE OBJECTIVE REAL-PARAMETER NUMERICAL OPTIMIZATION
# DEFINITIONS OF THE BENCHMARK FUNCTIONS BY JULIA
# REF: https://github.com/P-N-Suganthan/CEC2017-BoundContrained/blob/master/Definitions%20of%20%20CEC2017%20benchmark%20suite%20final%20version%20updated.pdf

# main caller of cec functions
# functionNumber is the number of the cec function which can be from 1 to 29
# x is the array of variables which will be used as parameters of the function
function cecFunction(functionNumber::Int64, x::Array{Float64, 1})::Float64

    if !(1 .< functionNumber .< 30)
        @error "functionNumber is not valid." functionNumber
        return false
    end

    theFunction = Symbol("cec",functionNumber)
    theFunction = @eval $theFunction
    return theFunction(x)

end

# BENT CIGAR FUNCTION
function cec1(x::Array{Float64, 1})::Float64
    
    term1 = x[1]^2
    term2 = 0
    for i in 2:length(x)
        term2 = term2 + x[i]^2
    end
    term2 = term2 * 1000000

    return term1 + term2

end

# ZAKHAROV FUNCTION
function cec2(x::Array{Float64, 1})::Float64

    return sum(x.^2) + sum(0.5 * x)^2 + sum(0.5 * x)^4

end

# ROSENBROCK'S FUNCTION
function cec3(x::Array{Float64, 1})::Float64

    term = 0
    for i in 1:length(x)-1
        term = term + 100 * (x[i]^2 - x[i+1])^2 + (x[i] - 1)^2
    end

    return term

end

# RASTRIGIN'S FUNCTION
function cec4(x::Array{Float64, 1})::Float64

    term = 0
    for i in 1:length(x)
        term = term + x[i]^2 - 10 * cos(2 * pi * x[i]) + 10
    end

    return term

end

# EXPANDED SCHAFFER'S F6 FUNCTION
function cec5(x::Array{Float64, 1})::Float64

    push!(x, x[1])

    term = 0
    for i in 1:length(x)-1    
    
        term = term + cec5g(x[i], x[i+1])

    end
    
    return term

end

function cec5g(x::Float64, y::Float64)::Float64

    term1 = sin(sqrt(x^2 + y^2))^2 - 0.5
    term2 = (1 + 0.001*(x^2 + y^2))^2

    return 0.5 + (term1 / term2)

end

# LUNACEK BI-RASTRIGIN FUNCTION
function cec6(x::Array{Float64, 1})::Float64

    d = 1
    L = length(x)
    D = L
    s = abs(1-(1/(2 * sqrt(D + 20) - 8.2)))
    m0 = 2.5
    m1 = -1 * sqrt((m0^2-d)/s)
    
    term1 = min(sum((x.-m1).^2) , d * D + s * sum((x.-m1).^2))

    term2 = 0
    for i in 1:L

        term2 = term2 + 1 - cos(2 * pi * (x[i] - m0))

    end

    return term1 + 10 * term2

end

# NON-CONTINUOUS ROTATED RASTRIGIN'S FUNCTION
function cec7(x::Array{Float64, 1})::Float64



end

# LEVY FUNCTION
function cec8(x::Array{Float64, 1})::Float64

    L = length(x)
    term = sin(pi * cec8w(x[1]))^2
    D = convert(Float64, L)
    
    for i in 1:(L - 1)

        term = term + 
        (cec8w(x[i]) - 1)^2 * (1 + 10 * sin(pi * cec8w(x[i]) + 1)^2) + 
        (cec8w(D) - 1)^2 * (1 + sin(2 * pi * cec8w(D))^2)

    end

    return term

end

function cec8w(x::Float64)::Float64

    return 1 + ((x - 1)/4)

end

# MODIFIED SCHWEFEL'S FUNCTION
function cec9(x::Array{Float64, 1})::Float64

    L = length(x)
    d = convert(Float64, L)
    term = 0

    for i in 1:L

        z = x[i] + 420.9687462275036
        term = term + cec9g(z, d)

    end

    return 418.9829 * d - term

end

function cec9g(z::Float64, d::Float64)::Float64

    if abs(z) <= 500

        return z * sin(abs(z)^(1/2))

    elseif z > 500

        return (500 - mod(z, 500)) * sin(sqrt(abs(500 - mod(z, 500)))) - ((z - 500)^2/10000 * d)

    else

        return (mod(abs(z),500)-500) * sin(sqrt(abs(mod(abs(z,500))-500))) - ((z + 500)^2/10000 * d)

    end

end

# HIGH CONDITIONED ELLIPTIC FUNCTION
function cec10(x::Array{Float64, 1})::Float64

    L = length(x)
    d = convert(Float64, L)
    term = 0
    for i in 1:L

        term = term + (10^6)^((i - 1)/(d - 1)) * x[i]^2

    end

    return term

end

# DISCUS FUNCTION
function cec11(x::Array{Float64, 1})::Float64

    term = 10^6 * x[1]^2

    for i in 2:length(x)

        term = term + x[i]^2

    end

    return term

end

# ACKLEY'S FUNCTION
function cec12(x::Array{Float64, 1})::Float64

    L = length(x)
    d = convert(Float64, L)
    term1 = 0
    for i in 1:L
        term1 = term1 + x[i]^2
    end

    term2 = 0
    for i in 1:L
        term2 = term2 + cos(2 * pi * x[i])
    end

    return -20 * exp(-0.2 * sqrt((1/d) * term1)) - exp((1/d)*term2) + 20 + Base.MathConstants.e

end

# WEIERSTRASS FUNCTION
function cec13(x::Array{Float64, 1})::Float64

    L = length(x)
    a = 0.5
    b = 3
    kmax = 20
    d = convert(Float64, L)

    term3 = 0
    for i in 1:L

        term1 = 0
        for k in 0:kmax

            term1 = term1 + a^k * cos(2 * pi * b^k * (x[i] + 0.5))

        end

        term2 = 0
        for k in 0:kmax

            term2 = term2 + a^k * cos(2 * pi * b^k * 0.5)

        end

        term3 = term3 + term1 - d * term2

    end

    return term3

end

# GRIEWANK'S FUNCTION
function cec14(x::Array{Float64, 1})::Float64

    L = length(x)

    term1 = 0
    for i in 1:L

        term1 = term1 + (x[i]^2 / 4000)

    end

    term2 = 1
    for i in 1:L

        term2 = term2 * cos(x[i]/sqrt(i))

    end

    return term1 - term2 + 1

end

# KATSUURA FUNCTION
function cec15(x::Array{Float64, 1})::Float64

    L = length(x)
    d = convert(Float64, L)
    term1 = 1
    for i in 1:L

        term2 = 0
        for j in 1:32

            term2 = term2 + (abs((2^j)*x[i] - round((2^j) * x[i]))/(2^j))

        end

        term1 = term1 * (1 + i * term2)^(10/(d^12))

    end

    return (10/(d^2)) * term1 - (10/(d^2))

end

# HAPPYCAT FUNCTION
function cec16(x::Array{Float64, 1})::Float64

    d = convert(Float64, length(x))
    L = length(x)

    term1 = 0
    for i in 1:L
        term1 = term1 + x[i]
    end

    term2 = 0
    for i in 1:L
        term2 = term2 + x[i]^2
    end

    return abs(term2 - d)^(1/4) + (0.5 * term2 + term1) / d + 0.5

end

# HGBAT FUNCTION
function cec17(x::Array{Float64, 1})::Float64

    L = length(x)
    d = convert(Float64, L)
    term1 = 0
    for i in 1:L
        term1 = term1 + x[i]
    end

    term2 = 0
    for i in 1:L
        term2 = term2 + x[i]^2
    end

    return abs(term2^2 - term1^2)^(1/2) + (0.5 * term2 + term1)/d + 0.5

end

# EXPANDED GRIEWANK'S PLUS ROSENBROCK'S FUNCTION
function cec18(x::Array{Float64, 1})::Float64

    # f7: cec6
    # f4: cec3

    L = length(x)
    term = 0
    for i in 1:(L - 1)

        term = term + cec6([cec3([x[i], x[i+1]])])

    end

    return term + cec6([cec3([x[L], x[1]])])

end

# SCHAFFER'S F7 FUNCTION
function cec19(x::Array{Float64, 1})::Float64

    L = length(x)
    term = 0
    for i in 1:(L-1)

        term = term + sqrt(cec19s(x[i], x[i+1])) * sin(50 * cec19s(x[i], x[i+1])^0.2)

    end

    return ((1.0 / (L - 1.0)) * term)^2

end

function cec19s(x::Float64, y::Float64)::Float64

    return sqrt(x^2 + y^2)

end