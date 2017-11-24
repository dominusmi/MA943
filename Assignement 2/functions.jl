# Recurrent function which trasverses a linked list and prints the data of each node
# @param:  list::Nullable{LList}
# @return: none
function list_trasverser(list::Nullable{LList})
    if isnull( list )
        # Check whether at leaf node
        return true;
    else
        # If not, show data
        node = get(list)
        key = node.data.key
        value = node.data.value
        println("$(key): $(value)")
        list_trasverser(node.next)
    end
end

# 'Linearly' searches a linked list for the node with a key value 'k'
# @param
# list: linked list to search
# k: key to search for
# @return
# ::Pair: Pairnode with key 'k' if found, othw 'false'
function search_list(list::Nullable{LList}, k::Int)
    if isnull( get(list).next )
        # Check whether at leaf node
        return Nullable{KVPair};
    else
        # If not, check current node key, and otherwise pass on
        node = get(list)
        if( node.data.key == k )
            return node.data
        else
            return search_list(node.next, k)
        end
    end
end

# Generates a linked list of random number of size 'size', with increasing unique keys in 1:'size'
# @param
# size: size of the linked list
# @return
# list::Nullable{LList}
function generate_llist(size::Int)
    arr = Array(KVPair,0)
    for i in 1:size
        tmp = KVPair(i,rand())
        push!(arr, tmp)
    end

    return buildLList(arr)
end

# Generates a linked list containing the cumulative sum
# @param:
# size: sisze of list to create
# rng: random generator
# @return
# Nullable{LList}: The generated list
# Array{floats}: the array containing the cumulative sum at every node
function generate_cumulative_llist(size::Int, rng)
    # Set number of intervals
    n=size

    # Generate a sample of n uniform random numbers in the interval [0,1]
    X = rand(rng, n);

    # Generate array of key-values
    pairs = Array(KVPair,0)
    cumsum = zeros(n+1)
    for i in 1:n
        cumsum[i+1] = cumsum[i] + X[i]
        push!(pairs, KVPair(i, cumsum[i+1]))

    end

    return buildLList(pairs), cumsum
end

# 'Linearly' searches linked list for the node representing the interval containing 'x'
# @params
# list: list to search
# x: float whose interval is being searched
# @return
# Pair: pair of key-value which contains x
function intervalmembership(list::Nullable{LList}, x::Float64)
    if isnull(list)
        return -1
    else
        node = get(list)
        if node.data.value > x
            return node.data
        else
            return intervalmembership_brute(node.next,x)
        end
    end
end

# Searches the Fenwick tree for the node representing the interval containing 'x'
# @param
# FT: Tree to search
# x: value to search
# @return
# Pair: pair object containing the value x
function intervalmembership(FT::Nullable{FTree}, x::Float64)
    if isnull(FT)
        # Really shouldn't happen ...
        return -1
    else
        # Normal case
        node = get(FT)
        if isnull( node.right ) && isnull( node.left )
            # If at leaf, return current
            return node.data
        else
            # Else check which tree to return
            y = get(node.left).data.value
            if x <= y
                return intervalmembership(node.left, x)
            else
                return intervalmembership(node.right, x-y)
            end
        end
    end
end

# Solves diffusion numerically using linear search over linked lists
# @param
# N: size of array
# @return
# X: interval
# P1: theoretical solution
# P: numerical solution

function brute_diffusion(N,T)

    L  = 10.0
    Nx = 201
    dx = 2.0*L/(Nx-1)
    X  = dx.*(-(Nx-1)/2:(Nx-1)/2)
    Y  = zeros(Int64,N)
    D  = 1.0
    t  = 0.0

    Drates = randexp(rng, N)
    Drates = [Drates;Drates]
    rates = Drates/2 * 1/(dx^2)

    rate_intervals = generate_cumulative_from_pairs( KVPair.( collect(1:2*N), rates ) )

    Σ = sum(rates)
    totalRate = Σ
    dt = 1.0/totalRate
    T=1.0

    # This is the main loop
    while t < T
        # Pick an event
        rate_selected = rand()*Σ
        k = intervalmembership(rate_intervals,rate_selected).key
        if k<=N
            hop = 1
            particleId = k
        else
            hop = -1
            particleId=k-N
        end
        Y[particleId]+=hop
        t+=dt
    end

    # Calculate the estimated density of particles
    P = zeros(Float64,length(X))
    for i in 1:length(Y)
        P[Y[i]+Int64((Nx-1)/2)+1]+=1/(N * dx)
    end

    # Calculate the theoretical density and compare
    function normal(x, D, t)
        return exp( -sqrt(2 / (D*t) ) * abs(x) ) / sqrt( 2*D*t )
    end
    P1 = normal.(X, D, T)

    return X, P1, P
end



# Same as previous, except using Fenwick Tree search

function fenwick_diffusion(N,T)
    L  = 10.0
    Nx = 201
    dx = 2.0*L/(Nx-1)
    X  = dx.*(-(Nx-1)/2:(Nx-1)/2)
    Y  = zeros(Int64,N)
    D  = 1.0
    t  = 0.0

    # Generates rates related variables
    Drates = randexp(N)
    Drates = [Drates;Drates]
    rates = ( Drates/2 ) * ( 1/(dx^2) )
    # σ = cumsum(rates)
    # Σ = σ[end]
    Σ = sum(rates)

    rates_pairs = KVPair.( collect( [1:N;N+1:2*N] ), rates )
    rate_tree = Nullable{FTree}(FTree(KVPair(0,0.0)))
    rate_tree = buildFTree(rate_tree, rates_pairs);


    totalRate = Σ
    dt = 1.0/totalRate
    T=1.0

    # This is the main loop
    while t < T
        # Pick an event
        rate_selected = Float64(rand(Uniform(0, Σ)))
        k = intervalmembership(rate_tree, rate_selected)
        k = k.key
        if k<=N
            hop = 1
            particleId = k
        else
            hop = -1
            particleId=k-N
        end
        Y[particleId]+=hop
        t+=dt
    end

    # Calculate the estimated density of particles
    P = zeros(Float64,length(X))
    for i in 1:length(Y)
        P[Y[i]+Int64((Nx-1)/2)+1]+=1/(N * dx)
    end

    # # Calculate the theoretical density and compare
    function normal(x, D, t)
        return exp( -sqrt(2 / (D*t) ) * abs(x) ) / sqrt( 2*D*t )
    end

    P1 = normal.(X, D, T)

    return X,P1,P
end

# Generates intervals and cumulative sums -- only used for testing
function generate_intervals(size, rng)

    n=size
    X = rand(rng, n)

    # Now calculate the array of partial sums
    values = Array{KVPair}(n)
    cumsum = zeros(n+1)

    for i in 1:n
        values[i] = KVPair(i,X[i])
        cumsum[i+1] = cumsum[i]+X[i]
    end

    return values,cumsum,rng
end

# Generates cumulative linked list from kvpais - only used for testing
function generate_cumulative_from_pairs(pairs::Array{KVPair})

    size = length(pairs)
    cumulatives = Array{KVPair}(size)

    last = 0.0

    for i in 1:size
        key = pairs[i].key
        value = pairs[i].value
        cumulatives[i] = KVPair( key, value + last )
        last += value
    end

    return buildLList(cumulatives)
end
